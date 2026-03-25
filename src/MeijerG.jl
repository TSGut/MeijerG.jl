module MeijerG

using HypergeometricFunctions: pFq
using SpecialFunctions: gamma, besselk

export meijerg, meijerg_reduce

include("core.jl")
include("reduce.jl")

end
