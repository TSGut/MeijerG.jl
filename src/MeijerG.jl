module MeijerG

using HypergeometricFunctions: pFq
using SpecialFunctions: gamma, besselk, loggamma, logabsgamma

export meijerg, meijerg_slater

include("slater.jl")
include("perturb.jl")
include("reduce.jl")

end
