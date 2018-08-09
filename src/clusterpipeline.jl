# Consensus and fine splitting functions

"""
    denoise(seqs; <keyword arguments> )

                    rough_radius::Float64 = 0.01,
                    fine_radius::Float64 = 1.0,
                    reassign::Bool = true,
                    clust_multiplier::Float64 = 1.0,
                    min_devs::Int64 = 6,
                    initial_devs::Int64 = 20,
                    user_min_clust::Int64 = 5,
                    k::Int64 = 6,
                    distfunc = euclidean,
                    center = mean,
                    verbose::Int64 = 1,
                    cycle_lim = 20,
                    triangle::Bool = false,
                    degap_param::Bool = false,
                    Polish::Bool = true,
                    p_value::Float64 = 0.01


#Arguments
- `seqs`: array containing the sequences to be operated on
- `rough_radius::Float64`: distance between SNPs used in fine clustering
- `fine_radius::Float64`: distance between SNPs used in rough clustering
- `reassign::Bool = true`: indicates whether or not to reassign the rough singletons
back to consensus sequences
- `clust_multiplier::Float64 = 1.0`: minimum non-error cluster size
- `min_devs::Int64 = 6`: number of high standard deviation kmer to fine cluster on
- `initial_devs::Int64 = 20`: initial amount of high standard deviation kmer to test
- `user_min_clust::Int64 = 5`: minimum size of a cluster, competes with clust_multiplier
- `k::Int64 = 6`: length of kmer sequences calculated
- `distfunc = euclidean`: type of distance function calculated on the kmer sequences
- `center = mean`: sets the center of a centroid to the mean of the centroid
- `verbose::Int64 = 1`: when true, the status of the program will be printed to 'stdout'
- `cycle_lim = 20`: maximum number of cycles allowed in clustering call
- `triangle::Bool = false`: controls whether or not the triangle inequality is used in calculations
- `degap_param::Bool = false`: controls whether degap is used in consensus calculations
- `Polish::Bool = false`: if true, draws a consensus again from singletons after reassignment
- `p_value::Float64 = 0.01`: statistical threshold for splitting decisions

Performs clustering and sequence consensus.

Order of process:
- Rough clustering
- Fine clustering
- Consensus of clusters

Recommended `radius` for rough clustering: 0.01
Recommended `radius` for fine clustering: 1.0

"""
function denoise(seqs; rough_radius::Float64=0.01, fine_radius::Float64=1.0, reassign::Bool=true,
                clust_multiplier::Float64=1.0, min_devs::Int64=6, initial_devs::Int64=20, user_min_clust::Int64=5,
                k::Int64=6, distfunc=euclidean, center=mean, verbose::Int64=1, cycle_lim=3, triangle::Bool=false,
                degap_param::Bool=false, Polish::Bool=true, p_value::Float64=0.01)

    #transforms input sequences into vectors of kmer frequencies
    kmer_vecs = [kmer_count(seq, k) for seq in seqs];

    fine_indices =  kmer_split_clustering(kmer_vecs, seqs, rough_radius=rough_radius, fine_radius=fine_radius,
                        reassign=reassign, clust_multiplier=clust_multiplier, min_devs=min_devs, initial_devs=initial_devs,
                        user_min_clust=user_min_clust, k=k, distfunc=distfunc, center=center, verbose=verbose,
                        cycle_lim=cycle_lim, triangle=triangle, p_value=p_value)

    return kmer_split_consensus(kmer_vecs, seqs, fine_indices, verbose=verbose, distfunc=distfunc, reassign=reassign,
                            k=k, degap_param=degap_param, user_min_clust=user_min_clust, Polish=Polish)
end

"""
    kmer_split_clustering(kmer_vecs, seqs; rough_radius::Float64 = 0.01, fine_radius::Float64 = 1.0,  reassign::Bool = true,
                        clust_multiplier::Float64 = 1.0, min_devs::Int64 = 6, initial_devs::Int64 = 20,
                        user_min_clust::Int64 = 5, k::Int64 = 6, distfunc = euclidean, center = mean,
                        verbose::Int64 = 1, cycle_lim = 20, triangle::Bool = false, p_value::Float64 = 0.01)

Returns the fine clustering result from the given sequences.

"""
function kmer_split_clustering(kmer_vecs, seqs; rough_radius::Float64=0.01, fine_radius::Float64=1.0, reassign::Bool=true,
                        clust_multiplier::Float64=1.0, min_devs::Int64=6, initial_devs::Int64=20, user_min_clust::Int64=5,
                        k::Int64=6, distfunc=euclidean, center=mean, verbose::Int64=1, cycle_lim=20, triangle::Bool=true,
                        p_value::Float64=0.01)


    #clusters the indices of input sequences into groups originating from the same templates
    return fine_clustering(kmer_vecs, rough_radius=rough_radius, fine_radius=fine_radius, clust_multiplier=clust_multiplier,
                           min_devs=min_devs, initial_devs=initial_devs, user_min_clust=user_min_clust, distfunc=distfunc,
                           center=center, verbose=verbose, cycle_lim=cycle_lim, triangle=triangle, p_value=p_value)
end

"""
    kmer_split_consensus(kmer_vecs, seqs, fine_indices; verbose::Int64=1, distfunc=euclidean, reassign::Bool=true,
                       k::Int64=6, degap_param::Bool=false, user_min_clust::Int64=5, Polish::Bool=true)

Returns the consensus of the sequences after clustering.

"""
function kmer_split_consensus(kmer_vecs, seqs, fine_indices; verbose::Int64=1, distfunc=euclidean, reassign::Bool=true,
                            k::Int64=6, degap_param::Bool=false, user_min_clust::Int64=5, Polish::Bool=true)
    if verbose > 0
        println("Finding consensus sequences...")
    end
    #convert the clusters of indices into clusters of sequences
    seq_clusters = [seqs[fine_indices[ind]] for ind in 1:length(fine_indices)]

    #determine consensus sequences that spawned each of the sequence clusters
    consensus = [consensus_seq(clust, degap_param = degap_param) for clust in seq_clusters]

    #adds back singletons from a preliminary clustering step in fine_clustering in order to obtain more
    #accurate frequencies
    if reassign
        if verbose > 0
            println("Reassigning singleton inputs to consensus sequences to improve accuracy of frequencies")
        end
        ave_dist_arr = Array{Float64, 1}(length(seqs))
        consensus_vecs =  [kmer_count(seq, k) for seq in consensus];
        num_C = length(consensus_vecs)
        num_vecs = length(kmer_vecs)
        C_arr = consensus_vecs
        cluster_indices = [[] for i in 1:num_C]
        for i in 1:num_vecs
            curr_input = kmer_vecs[i]
            dists = [distfunc(C_arr[j], curr_input) for j in 1:num_C]
            min_index = indmin(dists)
            push!(cluster_indices[min_index], i)
        end
        sizes = [length(cluster) for cluster in cluster_indices]
        arr = find(user_min_clust.<sizes) #consensus seqs with zero inputs disappear!
        if verbose > 0
            println("After reassigning, there are $(length(arr)) templates")
        end
        if Polish
            cluster_indices=cluster_indices[arr]
            sizes=sizes[arr]
            seq_clusters = [seqs[cluster_indices[ind]] for ind in 1:length(cluster_indices)]
            consensus = [consensus_seq(clust, degap_param = degap_param) for clust in seq_clusters]
            return consensus, sizes, cluster_indices
        else
            return consensus[arr], sizes[arr], cluster_indices[arr]
        end
    else
        return consensus, [length(clust) for clust in fine_indices], fine_indices
    end
end

"""
    get_ave_dists(fine_indices, original_seqs, consensus_seqs; distfunc=euclidean, center=median, k::Int64=6)

Finds the average distance using the declared distance function
"""
function get_ave_dists(fine_indices, original_seqs, consensus_seqs; distfunc=euclidean, center=median, k::Int64=6)
    kmer_vecs = [kmer_count(seq, k) for seq in original_seqs];
    num_seqs = length(consensus_seqs)
    ave_dist_arr = Array{Float64, 1}(num_seqs)
    consensus_vecs = [kmer_count(seq, k) for seq in consensus_seqs];
    for ind in 1:num_seqs
        arr = consensus_vecs[ind]
        clust = fine_indices[ind]
        ave_dist_arr[ind] = center([distfunc(arr, kmer_vecs[index]) for index in clust])
    end
    return ave_dist_arr
end

"""
    denoise_pipeline(infile::String; <Keyword Arguments> )

                         prefix = "",
                         outfile = "",
                         ref_file = "",
                         error_filter = 0.01,
                         homopoly_filter = true,
                         minlength = 20,
                         maxlength = 1000000,
                         trim_seqs = true,
                         align_ref_frames = true,
                         write_bad_seqs = false,
                         write_aminoacids = true,
                         minsize = 5,
                         cycle_lim = 10,
                         merge_differing_homopolymers = false,
                         cluster_filedump_path = "",
                         post_filedump = false,
                         rough_radius::Float64=0.01,
                         fine_radius::Float64=1.0,
                         reassign::Bool=true,
                         clust_multiplier::Float64=1.0,
                         min_devs::Int64=6,
                         initial_devs::Int64=20,
                         user_min_clust::Int64=5,
                         k::Int64=6,
                         distfunc=euclidean,
                         center=mean,
                         verbose::Int64=1,
                         triangle::Bool=false,
                         degap_param::Bool=false,
                         Polish::Bool=true,
                         p_value::Float64=0.01)

# Arguments
- `prefix = ""`: prefix for names of sequences in output file (ie in case of multiple time points, etc.).
- `outfile = ""`: default is "[infile name without extension]-clusters.fasta".
- `ref_file = ""`: panel of reference seqs -- if provided, used to orient and trim nearest sequences in input file.
- `error_filter = 0.01`: if the input is a ".fastq" file, any sequence with a site with error greater than this value is removed before clustering.
- `homopoly_filter = true`: if true, remove sequences with erroneous run-on homopolymer regions (or trim the sequences if this region is on an end) via hidden markov model inference.
- `minlength = 20`: removes sequences shorter than this prior to clustering.
- `maxlength = 1000000`: removes sequences longer than this prior to clustering.
- `trim_seqs = true`: if a reference panel is provided (`ref_file`) and this is true, locally trims consensus seqs to nearest reference sequence
- `align_ref_frames = true`: corrects codon reading frames of inferred consensus sequences.
- `write_bad_seqs = false`: write sequences with uncorrectable reading frames to file
- `write_aminoacids = true`: write amino acid translations of consensus sequences to file
- `minsize = 5`: throw away clusters/consensuses with fewer than this many reads.
- `cycle_lim = 10` : TODO
- `merge_differing_homopolymers = false`: aggregates clusters with consensus seqs that only differ by a single gap in a homopolymer region.
- `cluster_filedump_path = ""`: if nonempty, writes original sequences to output files determined by resulting cluster in given directory. This string may also include a file name prefix.
- `post_filedump`: writes inferred consensus seqs to file after clustering, before aligning and trimming.
- `rough_radius::Float64`: distance between SNPs used in fine clustering
- `fine_radius::Float64`: distance between SNPs used in rough clustering
- `reassign::Bool = true`: indicates whether or not to reassign the rough singletons
back to consensus sequences
- `clust_multiplier::Float64 = 1.0`: minimum non-error cluster size
- `min_devs::Int64 = 6`: number of high standard deviation kmer to fine cluster on
- `initial_devs::Int64 = 20`: initial amount of high standard deviation kmer to test
- `user_min_clust::Int64 = 5`: minimum size of a cluster, competes with clust_multiplier
- `k::Int64 = 6`: length of kmer sequences calculated
- `distfunc = euclidean`: type of distance function calculated on the kmer sequences
- `center = mean`: sets the center of a centroid to the mean of the centroid
- `verbose::Int64 = 1`: when true, the status of the program will be printed to 'stdout'
- `cycle_lim = 20`: maximum number of cycles allowed in clustering call
- `triangle::Bool = false`: controls whether or not the triangle inequality is used in calculations
- `degap_param::Bool = false`: controls whether degap is used in consensus calculations
- `Polish::Bool = false`: if true, draws a consensus again from singletons after reassignment
- `p_value::Float64 = 0.01`: statistical threshold for splitting decisions


A clustering pipeline for sequence files.
Order of processes:
- Error filter (if .fastq)
- Length filter
- Homopolymer filter
- Length filter again
- Orient (either to reference panel or inferred coarse cluster centroid)
- Rough clustering
- Fine clustering
- Consensus of clusters


Recommended `radius` for rough clustering: 0.01
Recommended `radius` for fine clustering: 1.0

`infile` and `ref_file` must be either .fasta or .fastq.
If file is .fastq, then removes sequences with greater than `error_filter`
error.


"""
function denoise_pipeline(infile::String;
                         prefix = "", outfile = "", ref_file = "",
                         error_filter = 0.01,
                         homopoly_filter = true,
                         minlength = 20,
                         maxlength = 1000000,
                         trim_seqs = true,
                         align_ref_frames = true,
                         write_bad_seqs = false,
                         write_aminoacids = true,
                         minsize = 5,
                         cycle_lim = 10,
                         merge_differing_homopolymers = false,
                         cluster_filedump_path = "",
                         post_filedump = false,
                         rough_radius::Float64=0.01,
                         fine_radius::Float64=1.0,
                         reassign::Bool=true,
                         clust_multiplier::Float64=1.0,
                         min_devs::Int64=6,
                         initial_devs::Int64=20,
                         user_min_clust::Int64=5,
                         k::Int64=6,
                         distfunc=euclidean,
                         center=mean,
                         verbose::Int64=1,
                         triangle::Bool=false,
                         degap_param::Bool=false,
                         Polish::Bool=true,
                         p_value::Float64=0.01)
    seqs = []
    if infile[end-5:end] == ".fasta"
        seqs, names = read_fasta_with_names_in_other_order(infile, seqtype=String)
        phreds = nothing
    elseif infile[end-5:end] == ".fastq"
        seqs, phreds, names = read_fastq(infile, seqtype=String)
        println("read $(length(seqs)) seqs")
        # maybe put this somewhere else?
        seqs, phreds, names = quality_filter(seqs, phreds, names, errorRate=error_filter)
        println("kept $(length(seqs)) seqs")
    else
        error("Unsupported file type: $(infile[end-5:end])")
    end
    if verbose > 0
        println(length(seqs), " total sequences read from ", infile)
    end

    if verbose > 0
        println("Filtering...")
    end

    minlength = max(minlength, 1)
    @time @sync seqs, phreds, names = length_filter(Array{String, 1}(seqs), phreds, names, minlength, maxlength)
    # TODO: write to specific filtered file option
    if verbose > 0
        println(length(seqs), " total sequences after length filter")
    end

    if homopoly_filter
        @time @sync seqs, phreds, names = homopolymer_filter(Array{String, 1}(seqs), phreds, names)
        # TODO: write to specific filtered file option
        if verbose > 0
            println(length(seqs), " total sequences after homopolymer filter")
        end

        @time @sync seqs, phreds, names = length_filter(Array{String, 1}(seqs), phreds, names, minlength, maxlength)
        # TODO: write to specific filtered file option
        if verbose > 0
            println(length(seqs), " total sequences after length filter")
        end
    end

    if length(ref_file) > 0
        @time @sync seqs, phreds, names = orient_to_refs_file(Array{String, 1}(seqs), phreds, names, ref_file)
    else
        # orient to coarse cluster centroid
        @time @sync seqs, phreds, names = orient_strands(Array{String, 1}(seqs), phreds, names, get_coarse_centroid(seqs))
    end
    # TODO: write to specific filtered file option
    if verbose > 0
        println(length(seqs), " total sequences after orienting")
    end

    if verbose > 0
        println("Clustering and Consensusing...")
    end

    #clusts, sizes, fine_indices = denoise(seqs, triangle=true, cycle_lim=cycle_lim, user_min_clust=10)
    clusts, sizes, fine_indices = denoise(seqs, rough_radius=rough_radius, fine_radius=fine_radius, reassign=reassign,
                clust_multiplier=clust_multiplier, min_devs=min_devs, initial_devs=initial_devs, user_min_clust=user_min_clust,
                k=k, distfunc=distfunc, center=center, verbose=1, cycle_lim=cycle_lim, triangle=triangle,
                degap_param=degap_param, Polish=Polish, p_value=p_value)

    if length(clusts) == 0
        warn("No consensus sequences found; try lowering minimum cluster size")
    end

    if length(prefix) > 0
        prefix = prefix * "_"
    end
    if outfile == ""
        outfile = infile[1:end-6] * "-clusters.fasta"
    elseif outfile[end-5:end] != ".fasta"
        println("Output filetype must be fasta, writing to default filename")
        outfile = infile[1:end-6] * "-clusters.fasta"
    end

    if verbose > 0
        println("Trimming and aligning...")
    end
    if trim_seqs && length(ref_file) > 0
        if post_filedump
            names = [prefix * "seq" * string(i) * "_" * string(sizes[i]) for i in 1:length(sizes)]
            write_fasta(outfile[1:end-6] * "-untrimmed.fasta", clusts, names=names)
        end
        clusts = trim_to_refs(Array{String, 1}(clusts), ref_file)
        tmp = merge_consensuses([[clusts[i], sizes[i]] for i in 1:length(sizes)])
        clusts, sizes = [t[1] for t in tmp], [t[2] for t in tmp]
    end

    if align_ref_frames
        if post_filedump
            names = [prefix * "seq" * string(i) * "_" * string(sizes[i]) for i in 1:length(sizes)]
            write_fasta(outfile[1:end-6] * "-uncorrected.fasta", clusts, names=names)
        end
        goodclsts, badclsts = align_reference_frames([clusts, sizes])
        clusts, sizes = goodclsts
        badclusts, badsizes = badclsts
    end

    names = [prefix * "seq" * string(i) * "_" * string(sizes[i]) for i in 1:length(sizes)]
    write_fasta(outfile, clusts, names=names)
    if write_aminoacids
        write_fasta(outfile[1:end-6] * "-aminoacids.fasta", [translate_to_aa(s[1:end-end%3]) for s in clusts], names=names)
    end
    if align_ref_frames && write_bad_seqs
        badnames = ["seq" * string(i) * "_" * string(badsizes[i]) for i in 1:length(badsizes)]
        write_fasta(outfile[1:end-6] * "-badseqs.fasta", badclusts, names=names)
    end

    if verbose > 0
        numseqs = length(sizes) > 0 ? sum(sizes) : 0
        medsize = length(sizes) > 0 ? median(sizes) : 0
        println("\nConsensus sequences written to ", outfile)
        println("Total sequences included: ", numseqs)
        println("Median cluster size: ", medsize)
        println("Number of clusters/consensus seqs: ", length(sizes))
    end
end

simple_clusterpipeline = denoise_pipeline



#==
    OLD PIPELINE
==#


"""
    cluster_and_consensus(seqs::Array{String}, radii; <keyword arguments>)

Cluster sequences, compute the consensus of each, and locally refine consensus sequences to corresponding reads.
Returns consensus sequences and cluster sizes for each.
`radii` may be a single float64 or an array of float64.

- `k = 6`: size of kmers computed in kmer vector for approximate dirichlet process clustering.
- `distfunc = corrected_kmer_dist(k)'`: distance function used in clustering.
- `center = mean`: method of computing cluster centroids.
- `shift = 7`: window size use to locally refine reference sequences.
- `minsize::Int = 10`: minimum size of cluster to return.
- `maxsize::Int = 100`: maximum number of sequences used to refine consensus sequence of cluster.
- `cycle_lims = fill(20, length(radii))'`: number of iterations of clustering.
- `filedump = false`: write original sequences to files separated by cluster.
- `filedump_path = ""`: path to write clusters of seqs.
- `names = nothing`: names of sequences to write to file. If not given, names are "seq1", etc.
- `phreds = nothing`: phreds of sequences to write to file. If given, writes ".fastq", else writes ".fasta".
- `verbose::Bool = false`: progress updates printed to std_out.
"""
function cluster_and_consensus(seqs::Array{String}, radii;
                                    k = 6,
                                    distfunc = corrected_kmer_dist(k),
                                    center = mean,
                                    shift = 7,
                                    minsize::Int = 10,
                                    maxsize::Int = 100,
                                    cycle_lims = fill(20, length(radii)),
                                    filedump = false,
                                    filedump_path = "",
                                    names = nothing,
                                    phreds = nothing,
                                    verbose::Bool = false)

    if verbose
        println("Computing kmer vectors...")
    end
    kmer_vecs = [kmer_count(seq, k) for seq in seqs]
    if verbose
        println("Clustering...")
    end
    μs, sizes, cluster_indices, centroid_indices = dp_centers(kmer_vecs, radii,
                                                              distfunc=distfunc,
                                                              center=center,
                                                              cycle_lims=cycle_lims,
                                                              verbose=verbose)
    filt = sizes .>= minsize
    final_sizes = sizes[filt]
    if verbose
        println("Keeping $(sum(filt)) clusters")
    end
    filt_indices = cluster_indices[filt]
    consensus_seqs = Array{String}(length(filt_indices))
    for i in 1:length(filt_indices)
        if verbose
            println("Computing consensus for cluster $i with size $(sizes[filt][i])")
        end

        consensus_seqs[i] = refine_ref(seqs[centroid_indices[filt][i]],
                                       sample(seqs[filt_indices[i]],
                                              min(maxsize, length(filt_indices[i])),
                                              replace=false),
                                       shift=shift)
        if filedump
            # write fasta for each split cluster - includes < minsize clusters
            ss = seqs[filt_indices[i]]
            nms = names==nothing ? ["seq$(iii)" for iii in 1:length(ss)] : names[filt_indices[i]]
            pds = phreds==nothing ? nothing : phreds[filt_indices[i]]
            if pds != nothing
                write_fastq(filedump_path * "clust-$(i).fastq", ss, pds, names=nms)
            else
                write_fasta(filedump_path * "clust-$(i).fasta", ss, names=nms)
            end
        end
    end
    # TODO: merge identical consensuses? (cant use merge_consensuses)
    return consensus_seqs, final_sizes
end

"""
    fine_cluster_and_consensus(seqs::Array{String}, radii; <keyword arguments>).

                                    polylen = 3,
                                    minsize::Int = 10,
                                    maxsize::Int = 100,
                                    cycle_lims = fill(20, length(radii)),
                                    merge_differing_homopolymers = false,
                                    filedump = false,
                                    filedump_path = "",
                                    names = nothing,
                                    phreds = nothing,
                                    verbose::Bool = false)

Cluster sequences, compute the consensus of each, and split into sub-clusters on above-threshold-discrepancies.
Returns consensus sequences and cluster sizes for each sub-cluster.
`radii` may be a single float64 or an array of float64.

- `k = 6`: size of kmers computed in kmer vector for approximate dirichlet process clustering.
- `distfunc = corrected_kmer_dist(k)'`: distance function used in clustering.
- `center = mean`: method of computing cluster centroids.
- `shift = 7`: window size use to locally refine reference sequences.
- `polylen`: min length of homopolymer region to be determined as a homopolymer region (see `peak_split`).
- `minsize::Int = 10`: minimum size of cluster to return.
- `cycle_lims = fill(20, length(radii))'`: number of iterations of clustering.
- `filedump = false`: write original sequences to files separated by cluster.
- `filedump_path = ""`: path to write clusters of seqs.
- `merge_differing_homopolymers = false`: if inferred consensus seqs only differ by a single gap in a homopolymer region, merge their clusters. (Note: this occurs after filedump)
- `names = nothing`: names of sequences to write to file. If not given, names are "seq1", etc.
- `phreds = nothing`: phreds of sequences to write to file. If given, writes ".fastq", else writes ".fasta".
- `verbose::Bool = false`: progress updates printed to std_out.
"""
function fine_cluster_and_consensus(seqs::Array{String}, radii;
                                    k = 6,
                                    distfunc = corrected_kmer_dist(k),
                                    center = mean,
                                    shift = 7,
                                    polylen = 3,
                                    minsize::Int = 10,
                                    cycle_lims = fill(20, length(radii)),
                                    merge_differing_homopolymers = false,
                                    filedump = false,
                                    filedump_path = "",
                                    names = nothing,
                                    phreds = nothing,
                                    verbose::Bool = false)

    if verbose
        println("Computing kmer vectors...")
    end
    kmer_vecs = [kmer_count(seq, k) for seq in seqs]
    if verbose
        println("Clustering...")
    end
    μs, sizes, cluster_indices, centroid_indices = dp_centers(kmer_vecs, radii,
                                                              distfunc=distfunc,
                                                              center=center,
                                                              cycle_lims=cycle_lims,
                                                              verbose=verbose)
    filt = sizes .>= minsize
    final_sizes = sizes[filt]
    if verbose
        println("Keeping $(sum(filt)) clusters")
    end
    filt_indices = cluster_indices[filt]

    consensus_seqs = []
    for i in 1:length(filt_indices)
        if verbose
            println("Computing consensus for cluster $i with size $(sizes[filt][i])")
        end
        # tmp array with only this cluster
        tmp_consensuses = []
        reads = []
        nms = names==nothing ? ["seq$(iii)" for iii in 1:length(filt_indices[i])] : names[filt_indices[i]]
        pds = phreds==nothing ? nothing : phreds[filt_indices[i]]

        recursive_split(seqs[filt_indices[i]],
                        tmp_consensuses,
                        centroid=seqs[centroid_indices[filt][i]],
                        shift=shift,
                        minreads=minsize,
                        polylen=polylen,
                        k=k,
                        reads=reads,
                        phreds=pds,
                        names=nms)
        if filedump
            for j in 1:length(reads)
                # check if phreds exist
                if reads[j][2] != nothing
                    write_fastq(filedump_path * "-clust-$(i)-split-$(j).fastq", reads[j][1], reads[j][2], names=reads[j][3])
                else
                    write_fasta(filedump_path * "-clust-$(i)-split-$(j).fasta", reads[j][1], names=reads[j][3])
                end
            end
        end
        # add these new consensuses to total list
        append!(consensus_seqs, tmp_consensuses)
    end
    if merge_differing_homopolymers
        consensus_seqs = merge_diff_homopolymers(consensus_seqs, polylen=polylen)
    else
        consensus_seqs = merge_consensuses(consensus_seqs)
    end
    filter!(clust->clust[2] >= minsize, consensus_seqs)
    return [c[1] for c in consensus_seqs], [c[2] for c in consensus_seqs]
end



#-----Maybe move these to different files------


"""
    trim_to_refs(seqs::Array{String, 1}, refspath::String)

Trims each sequence in `seqs` to nearest reference in `refspath` file after local
alignment. Distance to references determined by amino acid similarity (kmer dot product).
`refspath` may be either .fasta or .fastq.
"""
function trim_to_refs(seqs::Array{String, 1}, refspath::String)
    if length(refspath) < 5
        return seqs
    end
    if refspath[end-5:end] == ".fasta"
        refs = read_fasta_with_names_in_other_order(refspath, seqtype=String)
    elseif refspath[end-5:end] == ".fastq"
        refs = read_fastq(refspath, seqtype=String)
    else
        error("Unsupported file type: $(infile[end-5:end])")
    end
    refs = refs[1]
    return trim_to_refs(seqs, refs)
end

"""
    trim_to_refs(seqs::Array{String, 1}, refs::Array{String, 1}; kmersize = 6)

Trims each sequence in `seqs` to nearest reference in `refs` array after local
alignment. Distance determined by amino acid similarity (kmer dot product).
"""
function trim_to_refs(seqs::Array{String, 1}, refs::Array{String, 1}; kmersize = 6)
    if length(refs) == 0
        println("Not trimming: no reference sequences")
        return seqs
    end
    ref_kmers = [sparse_aa_kmer_count(rf, kmersize) for rf in refs]
    trimmed = seqs[:]
    function trim(sq::String)
        best_match = 0
        best_match_i = 0
        vec = sparse_aa_kmer_count(sq, kmersize)
        # compare to every reference, update min dists
        matches = [dot(vec, rfk) for rfk in ref_kmers]
        best_match_i = indmax(matches)
        if best_match_i > 0
            ali = loc_kmer_seeded_align(refs[best_match_i], sq)
            return degap(ali[2])
        else
            # this probably happens because it does not share any kmers with any references
            # fix by shortening kmersize?
            # make flag to throw these away?
            return sq
        end
    end
    trimmed = map(trim, trimmed)
    return trimmed
end
