module RAD

	include("include_all.jl")

	export

	# consensus.jl

	fine_clustering,
    gen_arrangements,
    cluster_split,
    get_non_homopolymers,
    are_homopolymers,

	# clusterpipeline.jl

	denoise,
    kmer_split_clustering,
    kmer_split_consensus,
    get_ave_dists,
    clusterpipeline,
		

end
