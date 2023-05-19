module SPTNRTLClapeyronExt
using SPTNRTL
using Clapeyron

"""
    sptNRTL::aspenNRTL

    function sptNRTL(smiles::Vector{String} components = smiles;
    puremodel=BasicIdeal,
    pure_userlocations = String[],
    verbose=false)

## Input parameters
 - none (provided by SPTNRTL.jl)

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
NRTL (Non Random Two Fluid) activity model, obtained from SPT-NRTL database. it uses `Clapeyron.aspenNRTL` to store the parameters.
```
Gᴱ = nRT∑[xᵢ(∑τⱼᵢGⱼᵢxⱼ)/(∑Gⱼᵢxⱼ)]
Gᵢⱼ exp(-αᵢⱼτᵢⱼ)
αᵢⱼ = αᵢⱼ₀ + αᵢⱼ₁T
τᵢⱼ = tᵢⱼ₀ + tᵢⱼ₁/T + tᵢⱼ₂*ln(T) + tᵢⱼ₃*T
```

## References
1. Renon, H., & Prausnitz, J. M. (1968). Local compositions in thermodynamic excess functions for liquid mixtures. AIChE journal. American Institute of Chemical Engineers, 14(1), 135–144. [doi:10.1002/aic.690140124](https://doi.org/10.1002/aic.690140124)
2. Winter, B., Winter, C., Esper, T., Schilling, J., & Bardow, A. (2023). SPT-NRTL: A physics-guided machine learning model to predict thermodynamically consistent activity coefficients. Fluid Phase Equilibria, 568(113731), 113731. [doi:10.1016/j.fluid.2023.113731](https://doi.org/10.1016/j.fluid.2023.113731)
"""
function SPTNRTL.sptNRTL(smiles::Vector{String},components = smiles; puremodel=Clapeyron.BasicIdeal,
    pure_userlocations = String[],
    verbose=false)

    a0,a1,t0,t1,t2,t3 = SPTNRTL.spt_NRTL_params(smiles)
    a0 = PairParam("a0",components,a0)
    a1 = PairParam("a1",components,a1)
    t0 = PairParam("t0",components,t0)
    t1 = PairParam("t1",components,t1)
    t2 = PairParam("t2",components,t2)
    t3 = PairParam("t3",components,t3)
    _puremodel = Clapeyron.init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = Clapeyron.aspenNRTLParam(a0,a1,t0,t1,t2,t3)
    references = String["10.1002/aic.690140124","10.1016/j.fluid.2023.113731"]
    model = Clapeyron.aspenNRTL(components,packagedparams,_puremodel,references)
    return model
end

end #module