module SPTNRTL

using CSV
using CSV.Tables
using Preferences
using PrecompileTools
using Downloads

#result type
struct SPTResult
    vals::NTuple{2,String}
    a::NTuple{2,Float64}
    tij::NTuple{4,Float64}
    tji::NTuple{4,Float64}
end

# cache, to avoid repeatingly obtain values
const SPT_NRTL_CACHE = Dict{NTuple{2,String},SPTResult}()

#source of SPT_NRTL database. in the format given by the paper.
const DEFAULT_NRTL_DB = "https://raw.githubusercontent.com/ClapeyronThermo/spt-nrtl-db/main/v3"
const SPT_NRTL_DB = Ref{String}(DEFAULT_NRTL_DB)

const DOWNLADED_TABLES = Dict{String,String}()

#to check if the source is online or not.
const ONLINE = Ref{Bool}(true)


const PERCENT_ENCODING_REPLACE = (
    " " => "%20",
    "!" => "%21",
    "\"" => "%22",
    "#" => "%23",
    "\$" => "%24",
    "%" => "%25",
    "&" => "%26",
    "'" => "%27",
    "(" => "%28",
    ")" => "%29",
    "*" => "%2A",
    "+" => "%2B",
    "," => "%2C",
    "/" => "%2F",
    ":" => "%3A",
    ";" => "%3B",
    "=" => "%3D",
    "?" => "%3F",
    "@" => "%40",
    "[" => "%5B",
    "]" => "%5D",
)

function percent_encoding(i::AbstractString)
    return replace(i,PERCENT_ENCODING_REPLACE...)
end

function _build_path(path,smiles_i,online = true)
    s1 = string(first(smiles_i))
    l1 = string(length(smiles_i))
    if online
        return path * '/' * s1 * '/' * l1 * '/' * percent_encoding(smiles_i) * ".csv"
    else
        return joinpath(abspath(normpath(path)),s1,l1,smiles_i * ".csv")
    end
end

function SPTResult(table,smiles_pair)
    smiles_2_list = Tables.columns(table)[:SMILES1]
    si,s2 = smiles_pair
    i = findfirst(isequal(s2),smiles_2_list)
    row = Tables.rows(table)[i]
    a1,a2 = row[:a_1],row[:a_2]
    a = (a1,a2)
    tij1,tij2,tij3,tij4 = row[:t_12_1],row[:t_12_2],row[:t_12_3],row[:t_12_4]
    tij = (tij1,tij2,tij3,tij4)
    tji1,tji2,tji3,tji4 = row[:t_21_1],row[:t_21_2],row[:t_21_3],row[:t_21_4]
    tji = (tji1,tji2,tji3,tji4)
    return SPTResult(smiles_pair,a,tij,tji)
end

function set_result!(mats,vals::SPTResult,smiles_pair,i,j) 
    a1,a2,t1,t2,t3,t4 = mats
    inverse = vals.vals != smiles_pair
    if inverse
        #check that we are actually inverse
        smiles_pair_inv = (last(smiles_pair),first(smiles_pair))
        true_inverse = vals.vals == smiles_pair_inv
        if !true_inverse
            @error ("incompatible SMILES pairs passed and somehow got into a set_result! function.")
            return nothing
        end
    end

    a = vals.a
    if !inverse
        tij = vals.tij
        tji = vals.tji
    else
        tij = vals.tji
        tji = vals.tij 
    end
    a1ij,a2ij = a
    t1ij,t2ij,t3ij,t4ij = tij
    t1ji,t2ji,t3ji,t4ji = tji
    a1[i,j] = a1ij
    a1[j,i] = a1ij
    a2[i,j] = a2ij
    a2[j,i] = a2ij

    t1[i,j] = t1ij
    t2[i,j] = t2ij
    t3[i,j] = t3ij
    t4[i,j] = t4ij
   
    t1[j,i] = t1ji
    t2[j,i] = t2ji
    t3[j,i] = t3ji
    t4[j,i] = t4ji
    return nothing
end

function set_db_path!(path,online = false) 
    @set_preferences!("path" => path)
    @set_preferences!("online" => online)
    @info "$path will be the default on the next session."
    ONLINE[] = online
    SPT_NRTL_DB[] = path
end

function db_path()
    read_preferences()
    SPT_NRTL_DB[]
end

function db_online()
    read_preferences()
    ONLINE[]
end
#initialization.
const INIT = Ref{Bool}(false)

function read_preferences()
    if !INIT[]
        path = @load_preference("path", DEFAULT_NRTL_DB)
        online = @load_preference("online", true)
        SPT_NRTL_DB[] = path
        ONLINE[] = online
        INIT[] = true
    end
    return nothing
end

function spt_NRTL_params(smiles_list;from = db_path(),online = db_online())
    n = length(smiles_list)
    n == 1 && throw(error("At least 2 SMILES are required."))
    a1 = zeros(Float64,(n,n))
    a2 = zeros(Float64,(n,n))
    t1 = zeros(Float64,(n,n))
    t2 = zeros(Float64,(n,n))
    t3 = zeros(Float64,(n,n))
    t4 = zeros(Float64,(n,n))
    mats = (a1,a2,t1,t2,t3,t4)
    found_vals = false
    for i in 1:n
        for j in 1:(i-1)
            smiles_pair = (smiles_list[i],smiles_list[j])
            smiles_pair_inv = (smiles_list[j],smiles_list[i])
            #check if the value is in the cache DB
            found_vals = haskey(SPT_NRTL_CACHE,smiles_pair)
            found_vals && set_result!(mats,SPT_NRTL_CACHE[smiles_pair],smiles_pair,i,j)
            found_vals && continue
            #check for the inverse
            found_vals = haskey(SPT_NRTL_CACHE,smiles_pair_inv)
            found_vals && set_result!(mats,SPT_NRTL_CACHE[smiles_pair_inv],smiles_pair_inv,i,j)
            found_vals && continue
            #retrieve CSV from DB
            src_path = _build_path(from,smiles_pair[1],online)
            if online & !found_vals
                DOWNLADED_TABLES
                #download table, store downloaded path for the session
                path = get!(DOWNLADED_TABLES,src_path,Downloads.download(src_path))
            else
                path = src_path
            end
            table = CSV.File(path)
            res = SPTResult(table,smiles_pair)
            
            #set result into the matrices
            set_result!(mats,res,smiles_pair,i,j)

            #store the cached result
            SPT_NRTL_CACHE[smiles_pair] = res
            SPT_NRTL_CACHE[smiles_pair_inv] = res
        end
    end
    return a1,a2,t1,t2,t3,t4
end

#Clapeyron.jl support, loaded as an extension
function sptNRTL end

export spt_NRTL_params,sptNRTL

@setup_workload begin
    path = normpath(Base.pkgdir(SPTNRTL),"test","test_db","db")
        @compile_workload begin
        spt_NRTL_params(["test1","test2"],from = path,online = false)
        spt_NRTL_params(["test1","test2"],from = path,online = false)
    end
    empty!(DOWNLADED_TABLES)
    empty!(SPT_NRTL_CACHE)
end
end #module
