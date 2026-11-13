using Test
using SpecialFunctions
using MeijerG

include("test_utilities.jl")
include("test_confluentjets.jl");
include("test_api_and_symmetry.jl");
include("test_formatter.jl");
include("test_core.jl");
include("test_contour_comparison.jl");
include("test_confluence_clusters.jl");
include("test_logratio_jet.jl");
include("test_higher_order_confluence.jl");
include("test_extended_lanes.jl");
include("test_meijerg_jet.jl");
include("test_domain.jl");
include("test_reductions.jl");
include("test_slater.jl");
include("test_external_references.jl");
