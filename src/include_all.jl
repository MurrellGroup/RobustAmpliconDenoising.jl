# For exporting names in imported packages to scripts that use BeNGS
# maybe remove for final version
using Reexport
using Compat
#using DPMeansClustering
#@reexport using NextGenSeqUtils
using Distributions

include("consensus.jl")
include("clusterpipeline.jl")
