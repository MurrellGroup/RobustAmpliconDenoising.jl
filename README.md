## Synopsis

Approximate dirichlet amplicon denoiser (RAD) contains cluster and consensus tools.

## Installation
```julia
Pkg.clone("RAD")

```

## Run tests
```julia
Pkg.test("RAD")

```

## Set paths
```julia
using RAD
```

<a id='Files-1'></a>
# Files
**`consensus.jl`** &mdash; *File*
```julia
fine_clustering(kmer_vecs; rough_radius::Float64 = 0.01, fine_radius::Float64 = 1.0, clust_multiplier = 1.0,
                     min_devs = 6, initial_devs = 20, user_min_clust = 10, distfunc = euclidean, center = mean,
                     verbose::Int64 = 1, cycle_lim = 20, triangle::Bool = false, p_value::Float64 = 0.01)

Preforms fine clustering on the already rough clustered sequences.
```

```julia
gen_arrangements(k)

Returns an array containing all possible kmer sequences of length k

```

```julia
get_non_homopolymers(dev_inds::Array{Int64, 1}, nuc_inds; min_devs = 6)

Takes in array of kmer indices. Returns the highest standard deviation non-homopolymer indices
```

```julia
are_homopolymers(alignment::Tuple{String,String}; k = 6, polylen = 3)

Determines if the aligned kmers represent a homopolymer difference. If the
kmers represent a homopolymer difference, the function returns true.

```


**`clusterpipeline.jl`** &mdash; *File*
```julia
clusterpipeline(infile::String; <Keyword Arguments> )

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

```

```julia
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
```

```julia
kmer_split_clustering(kmer_vecs, seqs; rough_radius::Float64 = 0.01, fine_radius::Float64 = 1.0,  reassign::Bool = true,
                    clust_multiplier::Float64 = 1.0, min_devs::Int64 = 6, initial_devs::Int64 = 20,
                    user_min_clust::Int64 = 5, k::Int64 = 6, distfunc = euclidean, center = mean,
                    verbose::Int64 = 1, cycle_lim = 20, triangle::Bool = false, p_value::Float64 = 0.01)

Returns the fine clustering result from the given sequences.
```

```julia
kmer_split_consensus(kmer_vecs, seqs, fine_indices; verbose::Int64=1, distfunc=euclidean, reassign::Bool=true,
                   k::Int64=6, degap_param::Bool=false, user_min_clust::Int64=5, Polish::Bool=true)

Returns the consensus of the sequences after clustering.
```

```julia
get_ave_dists(fine_indices, original_seqs, consensus_seqs; distfunc=euclidean, center=median, k::Int64=6)

Finds the average distance using the declared distance function
```
