module MeijerG

using HypergeometricFunctions: pFq
import Gamma
using Gamma: gamma, gamma_lower, gamma_inc, loggamma, logabsgamma, digamma
using Bessels: besselj, besselk, besseli
using LinearAlgebra

export meijerg, meijerg_slater, meijerg_laurent, meijerg_jet, meijerg_contour, @formatter, meijerg_to_latex_and_math

include("jets/confluentjets.jl")
using .ConfluentJets

include("utilities/errorfunctions.jl")
include("residues/slater.jl")
include("utilities/gammaratio.jl")
include("utilities/doubledouble.jl")
include("jets/confluentjets_compat.jl")
include("residues/confluence_clusters.jl")
include("residues/confluent_slater.jl")
include("jets/logratio_jet.jl")
include("jets/logratio_jet_extended.jl")
include("residues/slater_extended.jl")
include("residues/higher_laurent.jl")
include("jets/higher_logratio_jet.jl")
include("jets/meijerg_jet.jl")
include("contour/contourintegral.jl")
include("reductions/orthogonal_polynomials.jl")
include("reductions/reduce.jl")
include("core.jl")
include("formatter.jl")

end
