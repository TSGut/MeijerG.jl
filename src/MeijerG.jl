module MeijerG

using HypergeometricFunctions: pFq
using ClassicalOrthogonalPolynomials: laguerrel, jacobip
using SpecialFunctions: gamma, besselj, besselk, besseli, erf, erfc, loggamma, logabsgamma

export meijerg, meijerg_slater, meijerg_contour, @formatter, meijerg_to_latex_and_math

include("slater.jl")
include("contourintegral.jl")
include("perturb.jl")
include("reduce.jl")
include("core.jl")
include("formatter.jl")

end
