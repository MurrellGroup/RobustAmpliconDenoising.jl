
# Consensus and fine splitting functions


"""
   fine_clustering(kmer_vecs; rough_radius::Float64 = 0.01, fine_radius::Float64 = 1.0, clust_multiplier = 1.0,
                        min_devs = 6, initial_devs = 20, user_min_clust = 10, distfunc = euclidean, center = mean,
                        verbose::Int64 = 1, cycle_lim = 20, triangle::Bool = false, p_value::Float64 = 0.01)

Preforms fine clustering on the already rough clustered sequences.

"""

function fine_clustering(kmer_vecs; rough_radius::Float64=0.01, fine_radius::Float64=1.0, clust_multiplier=1.0,
                        min_devs=6, initial_devs=20, user_min_clust=10, distfunc=euclidean, center=mean, verbose::Int64=1,
                        cycle_lim=20, triangle::Bool=false, p_value::Float64=0.01)

    if verbose > 0
        println("Your inputs have begun rough clustering...")
    end
    if verbose == 1 ||verbose == 2
        v = 2
    elseif verbose == 0
        v = 0
    end

    #"rough" clustering step to get the massively different sequences into different clusters
    μs, sizes, indices, centroid_inds = cluster(kmer_vecs, rough_radius, distfunc=corrected_kmer_dist,
                                        center=center,verbose=v, cycle_lim=cycle_lim, triangle=triangle)



    #only sequences greater than indicated size may form a cluster
    keep_inds = findall(user_min_clust .< sizes)
    indices = indices[keep_inds]

    fine_indices = Array{Array{Int64, 1}, 1}()
    #get array which converts kmer indices into kmer sequences
    nuc_inds = gen_arrangements(Int(log(4, length(kmer_vecs[1]))))

    if verbose > 0
        println("Your inputs have begun fine clustering...")
    end
    if verbose == 1 || verbose == 0
        in_verbose = 0
    elseif verbose == 2
        in_verbose = 1
    end
    count=1

    #call cluster_split on every cluster created by rough clustering
    for ind in 1:length(indices)
        if verbose > 0
            println("Splitting cluster number $(ind)")
        end
        temp =cluster_split(count, kmer_vecs[indices[ind]], 1, nuc_inds, rough_radius=rough_radius, fine_radius=fine_radius,
            clust_multiplier=clust_multiplier, min_devs=min_devs, initial_devs=initial_devs, user_min_clust=user_min_clust,
            distfunc=euclidean,center=center, verbose=in_verbose, cycle_lim=cycle_lim, triangle=triangle, p_value=p_value)
        for i in 1:length(temp)
            push!(fine_indices, indices[ind][temp[i]])
            f=open("/home/mchernys/Documents/Murrell_lab/RAD_project/interm/testnew_ifines"
            *string(count)*".txt", "a")
            #arr=[join(a) for a in fine_indices[end]]
            str=join(fine_indices[end])
            write(f, str)
            close(f)
            count+=1
        end

    end
    if verbose > 0
        println("Returning $(length(fine_indices)) templates")
    end
    return fine_indices
end


"""
     gen_arrangements(k)

Returns an array containing all possible kmer sequences of length k


"""
function gen_arrangements(k)
    alph = [["A"],["C"],["G"],["T"]]
    current = copy(alph)
    for i in 1:k-1
        current = vcat([[vcat(j,letter) for letter in alph] for j in current]...)
    end
    return [join(i) for i in current]
end
"""
	cluster_split(kmer_vecs, count, nuc_inds; rough_radius::Float64=0.01, fine_radius::Float64=1.0, clust_multiplier=1.0,
                        min_devs=6, initial_devs=20, user_min_clust=10, distfunc=euclidean,
                        center=mean, verbose::Int64=1, cycle_lim=20, triangle::Bool=false, p_value::Float64=0.01)

Determines if passed cluster is a 'real' cluster and recursively splits non-error clusters
"""
function cluster_split(ind, kmer_vecs, count, nuc_inds; rough_radius::Float64=0.01, fine_radius::Float64=1.0, clust_multiplier=1.0,
                        min_devs=6, initial_devs=20, user_min_clust=10, distfunc=euclidean,
                        center=mean, verbose::Int64=1, cycle_lim=20, triangle::Bool=false, p_value::Float64=0.01)
    if verbose == 2
        println("______________________________________________")
        println("Current level of cluster_split recursion is $(count)")
        count+=1
    end

    k = Int(log(4, length(kmer_vecs[1])))
    num_vecs = length(kmer_vecs)
    #stdevs=Array{Float64, 1}()
    """"
    for  kmer_ind in 1:4^k
        arr=Array{Float64, 1}()
        for seq_ind in 1:length(kmer_vecs)
            push!(arr, kmer_vecs[seq_ind][kmer_ind])
        end
        if ind==81
            f=open("/home/mchernys/Documents/Murrell_lab/RAD_project/interm2/81stnew_arrs.txt", "a")
            str=join(arr)
            write(f, str)
            close(f)
            println(arr)

        end
        push!(stdevs, std(arr))
    end
    """
    #arr of each kmer's stdev throughout the kmer vectors
    aple=[[kmer_vecs[seq_ind][kmer_ind] for seq_ind in 1:length(kmer_vecs)] for kmer_ind in 1:4^k]

    stdevs = [std([kmer_vecs[seq_ind][kmer_ind] for seq_ind in 1:length(kmer_vecs)]) for kmer_ind in 1:4^k]


    #arr of indices of high->low stdev kmer indices
    dev_inds = reverse(sortperm(stdevs))[1:initial_devs]#setting the number of kmers high results in a lack of splitting



    #eliminate kmer indices that represent homopolymer differences
    dev_inds = get_non_homopolymers(dev_inds, nuc_inds, min_devs=min_devs)

    if length(dev_inds) == 0
        return [[i for i in 1:length(kmer_vecs)]]
    end
    #kmer counts of high stdev kmer indices
    top_kmers =[kmer[dev_inds] for kmer in kmer_vecs]
    #cluster kmers by high stdev indices
    us,sizes,indices,C_inds =cluster(top_kmers, fine_radius,distfunc=distfunc,center=center,
                                            verbose=verbose,cycle_lim=cycle_lim, triangle=triangle)

    #replace 2nd largest with each size to determine whether or not to split on it
    #p_val_arr = [(1-cdf(Poisson(maximum(sizes)*(error_rate/4)), sizes[i])) *sum(kmer_vecs[1]) for i in 1:length(sizes)]

    #minimum non-error cluster size. Estimate of the sum of poisson distributions
    #min_clust = clust_multiplier*.25*2.5*sqrt(rough_radius*sum(kmer_vecs[1])*2*log(maximum(sizes)))
    #minimum cluster size is the largest number between user input and the minimum non-error size
    keep_inds = findall(user_min_clust .< sizes)
    #keep_inds = Array{Int64, 1}()
    #for ind in first_keep_inds


    #keep_inds = [ind for ind in first_keep_inds if (1-cdf(Poisson(maximum(sizes)*(rough_radius/4)), sizes[ind])) *sum(kmer_vecs[1]) < p_value]

    #if there are no clusters to keep, add everything to the largest one
    if length(keep_inds) == 0 || length(keep_inds) == 1
        return [[i for i in 1:length(kmer_vecs)]]
    else
        final_indices = Array{Array{Int64, 1}, 1}()
        keep_centers = [center(kmer_vecs[indices[keep]]) for keep in keep_inds]

        for keep in keep_inds
            push!(final_indices, indices[keep])
        end

        for i in 1:length(indices) #go through every cluster index
            if !(i in keep_inds) #if that cluster is not one of the kept clusters
                C_arr = indices[i] #get the array of seq indices corresponding to the cluster
                for index in C_arr #go through every seq indix in the non-kept cluster
                    push!(final_indices[findmin([distfunc(kmer_vecs[index], center) for center in keep_centers])[2]],
                            index)
                end
            end
        end
    end

    #if the input has gotten this far, it must be split. Call cluster_split on each current cluster
    for clust_ind in 1:length(final_indices)
        outer_clust = copy(final_indices[clust_ind])
        new_clusts = cluster_split(ind, kmer_vecs[outer_clust], count, nuc_inds, rough_radius=rough_radius, fine_radius=fine_radius,
                                    min_devs=min_devs, distfunc=distfunc,center=center,verbose=verbose,
                                    cycle_lim=cycle_lim, triangle=triangle, p_value=p_value)
        for new_clust in new_clusts
            push!(final_indices, outer_clust[new_clust])
        end
        final_indices = deleteat!(final_indices, clust_ind)

    end
    return final_indices
end

"""
    get_non_homopolymers(dev_inds::Array{Int64, 1}, nuc_inds; min_devs = 6)

Takes in array of kmer indices. Returns the highest standard deviation non-homopolymer indices

"""
function get_non_homopolymers(dev_inds::Array{Int64, 1}, nuc_inds; min_devs=6)
    k=log(4, length(nuc_inds))
    len = length(dev_inds)
    dev_seqs = nuc_inds[dev_inds]
    to_remove = Array{Int64, 1}()
    #go through every high deviation kmer index
    for outer_ind in 1:len
        #make sure to git every combination only once
        for inner_ind in outer_ind + 1:len
            #align corresponding sequences and determine if they represent homopolymer differences
            if are_homopolymers(nw_align(dev_seqs[outer_ind], dev_seqs[inner_ind], mismatch_cost=-100.0), k=k)
                push!(to_remove, outer_ind, inner_ind)
            end
        end
    end
    #difference between all indices and the homopolymer indices is the non homopolymer indices
    new_dev_inds = setdiff(dev_inds, dev_inds[unique(to_remove)])

    if length(new_dev_inds) < min_devs #if there aren't enough indices, empty array says "don't split!"
        return []
    else
        return new_dev_inds[1:min_devs] #pick the top non-homopolymer kmer indices
    end
end

"""
    are_homopolymers(alignment::Tuple{String,String}; k = 6, polylen = 3)

Determines if the aligned kmers represent a homopolymer difference. If the
kmers represent a homopolymer difference, the function returns true.

"""
function are_homopolymers(alignment::Tuple{String,String};k=6,polylen::Int64=3)

    al1 = alignment[1]
    al2 = alignment[2]
    len = length(alignment[1])
    if len != k + 1 #homopolymers only differ by one
        return false
    end
    if polylen > k - 1
        return false
    end
    ind1 = findfirst(x->x=='-', alignment[1])#findfirst(isequal('-'), alignment[1])
    ind2 = findfirst(x->x=='-', alignment[2])#findfirst(isequal('-'), alignment[2])
    diff = ind1-ind2
    #determine if inds are only on edges
    #        XXXXX- (ind1)
    # (ind2) -XXXXX
    if abs(diff) == k #both gaps are on edges
        str = al1[1:(1+polylen)] * al2[1:(1+polylen)]
        if (length(unique(al1[1:(1+polylen)] * al2[1:(1+polylen)])) == 2) || (length(unique(al1[(len-polylen):len] * al2[(len-polylen):len])) == 2)
            return true
        else
            return false
        end
    elseif ind1 == 1 || ind2 == 1 #one gap on left edge, one on non-edge
        if ind1 == 1
            temp = ind2
        else
            temp = ind1
        end
        if length(unique(al1[1:(1+polylen)] * al2[1:(1+polylen)])) == 2 #edge box is true. THIS MEANS NO HOMOPOLY!!!
            return false

        else #edge box is false. May be a homopoly.
            if temp + polylen <= len && length(unique(al1[temp:(temp+polylen)]*al2[temp:(temp+polylen)])) == 2
                return true
            elseif temp - polylen >= 1 && length(unique(al1[(temp-polylen):temp]*al2[(temp-polylen):temp])) == 2
                return true
            end
        end
        if ind1 == 1
            al = al2
        else
            al=al1
        end
        #if you're here, it means that the non-edge gap is inside, not bordering, the homopoly diff
        if polylen + 2 < len
            if abs(diff) == 1 && length(unique(al[1:polylen+2])) == 2
                return true
            elseif abs(diff) == 2 && length(unique(al[2:polylen+3])) == 2 && length(unique(al1[[ind1, ind2]]*al2[[ind1, ind2]])) == 2
                return true
            end
        end
    else #one of the indices is k+1
        if ind1 == len
            temp = ind2
        else
            temp = ind1
        end
        if length(unique(al1[(len-polylen):len] * al2[(len-polylen):len])) == 2 #edge box is true. THIS MEANS NO HOMOPOLY!!!
            return false
        else #edge box is false. May be a homopoly

            if temp + polylen <= len && length(unique(al1[temp:(temp+polylen)]*al2[temp:(temp+polylen)])) == 2

                return true
            elseif temp - polylen >= 1 && length(unique(al1[(temp-polylen):temp]*al2[(temp-polylen):temp])) == 2
                return true
            end
        end
        if ind1 == len
            al=al2
        else
            al=al1
        end
        if polylen + 2 < len
            if abs(diff) == 1 && length(unique(al[len-polylen-1:len])) == 2
                return true
            elseif abs(diff) == 2 && length(unique(al[len-polylen-2:len-1])) == 2 && length(unique(al1[[ind1, ind2]]*al2[[ind1, ind2]])) == 2
                return true
            end
        end
    end
    #if the alignment did not meet any of these conditions, it must be a nonhomopolymer diff
    return false
    #=
    arr = [alignment[1]... , alignment[2]... ]
    #max homopolymer length difference is 2
    if length(findin(arr, '-')) > 2
        return false
    end

    len = length(alignment[1])

    #check left->right
    for i in 1:len
        if i+polylen > len #no chance of sufficiently large polymer, break
            break;
        end
        if (alignment[1][i] == '-') #test if this is start of homopolymer
            @inbounds temp = alignment[2][(i+1):(i+polylen)]
            temp2 = alignment[1][i:i+polylen]

            if length(unique(temp)) == 1 && length(unique(temp2))== 2
                return true
            end
        elseif (alignment[2][i] == '-')
            @inbounds temp = alignment[1][(i+1):(i+polylen)]
            temp2 = alignment[2][i:i+polylen]
            if length(unique(temp)) == 1 && length(unique(temp2))== 2
                return true
            end
        end
    end
    #check right->left
    for i in len:-1:1
        if i - polylen < 1 #no chance of sufficiently large polymer, return
            return false
        end
        if (alignment[1][i] == '-')
            @inbounds temp = alignment[2][(i-polylen):(i-1)]
            temp2 = alignment[1][i-polylen:i]
            if length(unique(temp)) == 1 && length(unique(temp2)) == 2
                return true
            end
        elseif (alignment[2][i] == '-')
            @inbounds temp = alignment[1][(i-polylen):(i-1)]
            temp2 = alignment[2][i-polylen:i]
            if length(unique(temp)) == 1 && length(unique(temp2)) == 2
                return true
            end
        end
    end
    return false
    =#
end


#==
      OLD VERSION
==#

"""
    alignment_consensus(seqs)

Computes a consensus sequence for array of given sequences that have already been aligned.
Consensus is computed by site-wise mode among given sequences.
"""
function alignment_consensus(seqs)
    string_lengths = [length(i) for i in seqs]
    ncols = length(seqs[1])
    if !all([x == ncols for x in string_lengths])
        error("Aligned sequences must all be the same length.")
    end
    result = string([mode([seq[i] for seq in seqs]) for i in 1:ncols]...)
    return degap(result)
end

"""
    coords(ref, read)

Returns an array of degapped coordinates, such that
coords(ref, read)[i] gives you the position the aligned read/ref
that matches the i'th ungapped position in ref.
"""
function coords(ref, read)
    if length(ref) != length(read)
        error("Aligned strings are meant to be the same length.")
    end
    degappedRef = degap(ref)
    coordMap = Array{Int64, 1}(undef, length(degappedRef))
    count = 1
    for i in 1:length(degappedRef)
        while ref[count] == '-'
            count += 1
        end
        coordMap[i] = count
        count += 1
    end
    return coordMap
end

"""
    get_centroid(reads, k, distfunc = corrected_kmer_dist)

Computes a rough cluster centroid by returning the sequence in `reads`
that has the nearest distance in `k`-mer space to the mean `k`-mer vector.
"""
function get_centroid(reads, k, distfunc = corrected_kmer_dist)
    kmerVecs = [kmer_count(read, k) for read in reads]
    meanVec = mean(kmerVecs)
    dists = [distfunc(kmerVec, meanVec, k=k) for kmerVec in kmerVecs]
    return reads[findmin(dists)[2]]
end

"""
    consensus_seq(reads; thresh = 0.7, shift = 1, k = 6,
                       distfunc = corrected_kmer_dist)

Computes a cluster centroid by taking the nearest read to the mean `k`-mer
vector (see `get_centroid(reads, k, distfunc)`) and refining it (see `refine_ref(ref, reads)`).
`shift` determines the window size of comparison between sequences when refining rough centroid
locally (actual window size is `shift+1`).
"""
function consensus_seq(ind, reads; degap_param = true, thresh = 0.7, shift = 1, k = 6,
                       distfunc = corrected_kmer_dist)
    ref = get_centroid(reads, k, distfunc)
    return refine_ref(ref, reads, degap_param = degap_param, thresh = thresh, shift = shift)
end

#=function get_matches(candidate_ref, reads, shift; kmer_align = true)
    alignments = []
    if kmer_align
        alignments = [kmer_seeded_align(candidate_ref, i) for i in reads]
    else
        alignments = [nw_align(candidate_ref, i) for i in reads]
    end
    maps = [coords(i...) for i in alignments]
    matchContent = [[degap(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
    matches = [freq(matchContent[k], degap(candidate_ref[k:k+shift])) for k in 1:length(matchContent)]
    return alignments, maps, matches, matchContent
end=#
#warn("Multithreaded version active as sole get_matches")

"""
    get_matches(candidate_ref, reads, shift; kmer_align = true)

Determine alignments between each sequence of `reads` and `candidate_ref`, along with coordinate mappings
and local matches. Used internally in `refine_ref` and `peak_split`.
"""
function get_matches(candidate_ref, reads, shift; degap_param = true, kmer_align = true)
    alignments = []
    if kmer_align
        alignments = map(i -> kmer_seeded_align(candidate_ref, i), reads)
    else
        alignments = map(i -> nw_align(candidate_ref, i), reads)
    end

    maps = [coords(i...) for i in alignments]


    if (degap_param)
        matchContent = [[degap(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]

        matches = [freq(matchContent[k], degap(candidate_ref[k:k+shift])) for k in 1:length(matchContent)]

    else
       matchContent = [[(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]

       matches = [freq(matchContent[k], candidate_ref[k:k+shift]) for k in 1:length(matchContent)]


    end


    return alignments, maps, matches, matchContent
end

"""
    refine_ref(candidate_ref, reads;  thresh = 0.7, shift = 1)

Takes a candidate consensus sequence, `candidate_ref`, and corresponding array of `reads` and refines the candidate via majority
votes from reads at local windows of size `shift+1`. If after alignment the frequency of a local region of the candidate
is less than `thresh`, this part of the candidate is refined.

"""
function refine_ref(candidate_ref, reads; degap_param = true, thresh = 0.7, shift = 1)
    # need to add in code to handle deletions at the start and end of
    # the candidate!

    alignments, maps, matches, _ = get_matches(candidate_ref, reads, shift, degap_param = degap_param)

    boolVec = trues(length(candidate_ref))
#BitArray{1}(1+zeros(length(candidate_ref)))
    for i in 1:length(matches)
        if matches[i] < thresh
            for ind in i:i+shift
                boolVec[ind]=false
            end
        end
    end

    # need to explicitly handle the case where the reference is shorter than the reads.
    # refine front
    frontStrCol = [degap((alignments[j][2])[1:maps[j][1]-1]) for j in 1:length(maps)]
    frontNonEmpties = [length(s) != 0 for s in frontStrCol]
    frontConsensus = ""
    majority = length(alignments) / 2
    while length(findall(x->x!=0, frontNonEmpties)) > majority
        frontStrCol = frontStrCol[frontNonEmpties]
        frontConsensus = mode(frontStrCol)
        frontAli, frontMaps, _, _ = get_matches(frontConsensus, frontStrCol, 1, degap_param = degap_param, kmer_align=false)
        frontStrCol = [degap(frontAli[j][2][1:frontMaps[j][1]-1]) for j in 1:length(frontMaps)]
        frontNonEmpties = [length(s) != 0 for s in frontStrCol]
        majority = length(frontAli) / 2
    end
    frontStr = frontConsensus

    # refine end
    endStrCol = [degap((alignments[j][2])[maps[j][end]+1:end]) for j in 1:length(maps)]
    endNonEmpties = [length(s) != 0 for s in endStrCol]
    endConsensus = ""
    majority = length(alignments) / 2
    while length(findall(x->x!=0, endNonEmpties)) > majority
        endStrCol = endStrCol[endNonEmpties]
        endConsensus = mode(endStrCol)
        endAli, endMaps, _, _ = get_matches(endConsensus, endStrCol, 1, degap_param = degap_param, kmer_align=false)
        endStrCol = [degap(endAli[j][2][endMaps[j][end]+1:end]) for j in 1:length(endMaps)]
        endNonEmpties = [length(s) != 0 for s in endStrCol]
        majority = length(endAli) / 2
    end
    endStr = endConsensus

    # # if something breaks, CHECK THIS FIRST!
    if (length(boolVec)-sum(boolVec)) < 0.5
        return frontStr * candidate_ref * endStr
    end
    ranges = []
    i=1
    while i < length(boolVec)
        if !boolVec[i]
            start = i
            i += 1
            while !boolVec[i] && (i < length(boolVec))
                i += 1
            end
            push!(ranges,[start, i-1])
        end
        i += 1
    end

    newRef = candidate_ref[1:ranges[1][1]-1]
    for i in 1:length(ranges)-1
        low, high = ranges[i][1], ranges[i][2]
        insStr = mode([degap((alignments[j][2])[maps[j][low]:maps[j][high]])
                       for j in 1:length(maps)])
        newRef = newRef*insStr
        newRef = newRef*candidate_ref[ranges[i][2]+1:ranges[i+1][1]-1]
    end

    low, high = ranges[length(ranges)][1], ranges[length(ranges)][2]
    insStr = mode([degap((alignments[j][2])[maps[j][low]:maps[j][high]]) for j in 1:length(maps)])
    newRef = newRef*insStr
    newRef = newRef*candidate_ref[ranges[length(ranges)][2]+1:length(candidate_ref)]
    return frontStr*newRef*endStr
end

"""
    recursive_split(cluster, glom; centroid = "", shift = 7, minreads = 5, polylen = 3, k = 6, reads = nothing, phreds = nothing, names = nothing)

Recursively splits input `cluster` using `peak_split` and populates `glom` with the separated sub-clusters.
Splits clusters according to `centroid` (unless not given, in which case it is inferred).
`polylen` is lower bound for run-on subsequence to by considered a homopolymer region and not be split on.
Does not split on clusters of size less than `minReads`.
If `reads` is given as a string array, this will populate `reads` with tuples of `(cluster, phreds, names)` for
determining the actual clusters that are split, rather than just corresponding consensus sequences, where `phreds`
and `names` may be arrays of corresponding values, or `nothing`.
"""
function recursive_split(cluster::Array{String}, glom; centroid = "", shift = 7, minreads = 5, polylen = 3, k = 6, reads = nothing, phreds = nothing, names = nothing)
    consensus = ""
    if length(centroid) == 0
        consensus = consensus_seq(cluster, k=k)
    else
        consensus = refine_ref(centroid, cluster)
    end
    splitClus = peak_split(consensus, cluster, shift=shift, minReads=minreads, polylen=polylen)

    if length(splitClus) > 1
        for clust in splitClus
            split_phreds = phreds==nothing ? nothing : phreds[clust]
            split_names = names==nothing  ? nothing : names[clust]
            recursive_split(cluster[clust], glom, shift=shift, minreads=minreads, polylen=polylen, k=k, reads=reads, phreds=split_phreds, names=split_names)
        end
    else
        push!(glom, [consensus, length(cluster)])
        if reads != nothing
            push!(reads, (cluster, phreds, names))
        end
    end
end

"""
    peak_split(candidate_ref, reads; thresh = 0.7, shift = 3, minReads = 5, polylen = 3)

Creates and returns sub-clusters of `reads` split by variation from `candidateRef` if, locally, the most common read is less frequent than given `thresh`, default 0.7.
If the size of the input cluster `reads` is large, then `thresh` may be changed: thresh = max(1-10/length(reads), thresh).
Size of local window/site for comparison = `shift+1`.
Return values are arrays of boolean arrays, where each boolean array can be used to index into the input reads data to give the split clusters.
Determines  treshold for a cluster split, where if the candidate has local region with frequency less than `thresh` the cluster is split,
 to be `thresh = max(1-10/length(reads), thresh)` to account for larger cluster sizes.
`polylen` is lower bound for run-on subsequence to by considered a homopolymer region and not be split on.
Does not split on clusters of size less than `minReads`.
"""
function peak_split(candidate_ref, reads; thresh = 0.7, shift = 3, minReads = 5, polylen = 3)
    candidateRef = candidate_ref
    threshProp = thresh
    if length(reads) < minReads
        return [BitArray([true for r in reads])]
    end
    # Adjusting thresh is meant to account for the fact that 30% of a very large cluster is still very large,
    #  so we assume that, for large clusters, if we have 10 seqs, that is enough to split into a new cluster.
    # However, this will fail for very very large clusters that have 10 sequences that have the same mismatch base error.
    # Need to build in the expected number of mismatched base errors at any given base to control the false positive rate here.
    thresh = max(1-10/length(reads), thresh)
    numSteps = length(candidateRef) - shift
    #Need to add in code to handle deletions at the start and end of the candidate!
    alignments, maps, matches, matchContent = get_matches(candidateRef, reads, shift)
    boolVec = trues(length(candidateRef))
    for i in 1:length(matches)
        if (matches[i] < thresh)
            boolVec[i:i+shift] = false
        end
    end
    # TODO: make ranges statically sized to avoid push!()
    ranges = []
    i = 1
    while i <= numSteps
        if !boolVec[i]
            start = i
            while !boolVec[i] && i < numSteps
                i += 1
            end
            # Julia subranges are inclusive, => window size = (shift+1)
            if shift+1 < polylen
                start = max(1, start - (polylen - (shift+1)))
                i = min(numSteps, i + (polylen - (shift+1)))
            end
            push!(ranges, [start, i-1])
        end
        i += 1
    end

    # add 2 to lengths to account for front and end as well
    scoreList = 1 + zeros(length(ranges) + 2)
    scoreInds = zeros(UInt32, length(ranges) + 2)
    for i in 1:length(ranges)
        low, high = ranges[i][1], ranges[i][2]
        sortedFreqs = sorted_freqs([degap(alignments[j][2][maps[j][low]:maps[j][high]]) for j in 1:length(maps)])
        # The following might happen with discrepancies at ends of reads; this is taken care of later in this function
        if length(sortedFreqs) == 1
            continue
        end
        # This will not handle the case where the third most frequent variant is a true mutation and the second is noise, but this case only occurs or is important when:
        #  1) the second most frequent, noisy read is more common than reads of a mutation (ie in a very long homopolymer)
        #  2) the third variant, the true mutation, doesn't cause a split ANYWHERE else, so that the mutation is basically a single base change at a homopolymer region (as opposed to a noisy indel)
        # So it is very uncommon for this case to happen but it might need to be handled.
        ali = nw_align(sortedFreqs[1][2], sortedFreqs[2][2])
        if !differ_by_just_one_gap(ali, polylen=polylen)
            scoreList[i] = minimum(matches[ranges[i][1]:ranges[i][2]])
            scoreInds[i] = findmin(matches[ranges[i][1]:ranges[i][2]])[2] + ranges[i][1] - 1
        end
    end

    # take care of splits at overhanging beginnings of reads, where consensus is empty ("")
    # if any overhang is more common than threshold, potentially split on most common overhang
    frontStrCol = [degap(alignments[j][2][1:maps[j][1]-1]) for j in 1:length(maps)]
    frontOverhangs = [length(s) != 0 for s in frontStrCol]
    if freq(frontStrCol, "") < thresh
        sortedFreqs = sorted_freqs([degap(alignments[j][2][1:maps[j][polylen]])
                    for j in 1:length(maps) if frontOverhangs[j]])
        if length(sortedFreqs) == 1
            scoreList[end-1] = freq(frontStrCol, mode(frontStrCol[frontOverhangs]))
            scoreInds[end-1] = 1
        else
            ali = nw_align(sortedFreqs[1][2][1:min(polylen, end)], sortedFreqs[2][2][1:min(polylen, end)])
            if !differ_by_just_one_gap(ali, polylen=polylen)
                scoreList[end-1] = freq(frontStrCol, mode(frontStrCol[frontOverhangs]))
                scoreInds[end-1] = 1
            end
        end
    end
    # take care of splits at overhanging ends
    endStrCol = [degap(alignments[j][2][maps[j][end]+1:end]) for j in 1:length(maps)]
    endOverhangs = [length(s) != 0 for s in endStrCol]
    if freq(endStrCol, "") < thresh
        sortedFreqs = sorted_freqs([degap(alignments[j][2][maps[j][end-polylen]:end])
                    for j in 1:length(maps) if endOverhangs[j]])
        if length(sortedFreqs) == 1
            scoreList[end] = freq(endStrCol, mode(endStrCol[endOverhangs]))
            scoreInds[end] = 1
        else
            ali = nw_align(sortedFreqs[1][2][max(end-polylen, 0):end], sortedFreqs[2][2][max(end-polylen, 0):end])
            if !differ_by_just_one_gap(ali, polylen=polylen)
                scoreList[end] = freq(endStrCol, mode(endStrCol[endOverhangs]))
                scoreInds[end] = 1
            end
        end
    end
    # TODO: maybe add singleton clusters back to cluster with closest distance
    if minimum(scoreList) < thresh
        ind = findmin(scoreList)[2]
        splitpos = scoreInds[ind]
        # if in middle of sequence
        if ind < length(scoreInds) - 1
            unio = union(matchContent[splitpos])
            splitCluster = [matchContent[splitpos] .== s for s in unio]
            return splitCluster
        # else if on front
        elseif ind == length(scoreInds) - 1
            unio = union(frontStrCol)
            splitCluster = [frontStrCol .== s for s in unio]
            return splitCluster
        # else if on end
        else
            unio = union(endStrCol)
            splitCluster = [endStrCol .== s for s in unio]
            return splitCluster
        end
    else
        return [BitArray([true for r in reads])]
    end
end

"""
    differ_by_just_one_gap(alignment; polylen = 3)

`alignment` is an array containing two already aligned strings.
Returns true if first and second sequences in given array differ by
just one gap where that gap is in a homopolymer region.
A homopolymer region is determined to be a region of a single repeated
nucleotide of length at least `polylen`.
"""
function differ_by_just_one_gap(alignment; polylen = 3)
    gapDiffCount = 0
    otherDiffCount = 0
    gapInd = 0
    noGapSeq = 0
    for i in 1:length(alignment[1])
        if (alignment[1][i] == '-')
            gapDiffCount += 1
            gapInd = i
            noGapSeq = 2
        elseif (alignment[2][i] == '-')
            gapDiffCount += 1
            gapInd = i
            noGapSeq = 1
        elseif (alignment[1][i] != alignment[2][i])
            otherDiffCount += 1
        end
    end
    if otherDiffCount == 0 && gapDiffCount == 1
        # check if homopolymer region by matching regex to local string
        lo, hi = max(gapInd - (polylen-1), 1), min(gapInd + (polylen-1), length(alignment[1]))
        regex = Regex(string(alignment[noGapSeq][gapInd]) * "{" * string(polylen) * "}")
        return ismatch(regex, alignment[noGapSeq][lo:hi])
    else
        return false
    end
end

"""
    consensus_viz(candidate_ref, reads; thresh = 0.7, shift = 3,
                  intitle = "Consensus agreement.")

Creates a plot to visualize the agreement of a consensus sequence, `candidate_ref`,
with its cluster of `reads`.
Size of local window/site for comparison = `shift+1`.
I'm pretty sure `thresh` does nothing here.
"""
function consensus_viz(candidate_ref, reads;
                       thresh = 0.7, shift = 3, intitle = "Consensus agreement.")
    #Need to add in code to handle deletions at the start and end of the candidate!
    alignments, maps, matches, _ = get_matches(candidate_ref, reads, shift)
    frontStrCol = freq([degap((alignments[j][2])[1:maps[j][1]-1]) for j in 1:length(maps)], "")
    endStrCol = freq([degap((alignments[j][2])[maps[j][end]+1:end]) for j in 1:length(maps)], "")
    w, h = plt[:figaspect](0.1)
    figure(figsize=(w, h))
    plot(1:length(matches)+2, vcat([frontStrCol], matches,[endStrCol]))
    ax = gca()
    ax[:set_xlim]([-10, length(matches)+12])
    ax[:set_ylim]([0, 1])
    title(intitle)
    xlabel("Consensus Position")
    ylabel("Agreement")
end

"""
    disagreements(candidate_ref, reads; thresh = 0.7, shift = 3)

Prints local disagreements between `candidate_ref` and each sequence of `reads` after aligning.
Size of local window/site for comparison = `shift+1`.
A disagreement is a region where the candidate has a local region with frequency less than `thresh`.
"""
function disagreements(candidate_ref, reads; thresh = 0.7, shift = 3)
    # need to add in code to handle deletions at the start and end of
    # the candidate!
    alignments = [kmer_seeded_align(candidate_ref, i) for i in reads]
    maps = [coords(i...) for i in alignments]

    matches = [freq([(alignments[i][2])[maps[i][k]:maps[i][k+shift]]
                     for i in 1:length(maps)], candidate_ref[k:k+shift])
               for k in 1:length(candidate_ref)-shift]

    boolVec = trues(length(candidate_ref))
    for i in 1:length(matches)
        if matches[i] < thresh
            boolVec[i:i+shift] = boolVec[i:i+shift] & [false for k in 1:1+shift]
        end
    end

    ranges = []
    i = 1
    while i < length(boolVec)
        if !boolVec[i]
            start = i
            i += 1
            while !boolVec[i] && i < length(boolVec)
                i += 1
            end
            push!(ranges,[start, i-1])
        end
        i += 1
    end

    frontStrCol = [degap((alignments[j][2])[1:maps[j][1]-1]) for j in 1:length(maps)]
    if freq(frontStrCol, "") < thresh
        println("Before ref starts:")
        freq_dict_print(countmap(frontStrCol), thresh=1.1)
    end

    for i in 1:length(ranges)
        low, high = ranges[i][1], ranges[i][2]
        insStrCol = [degap((alignments[j][2])[maps[j][low]:maps[j][high]]) for j in 1:length(maps)]
        if freq(insStrCol, candidate_ref[low:high]) < thresh
            println("Ref pos ", low, ":", candidate_ref[low:high])
            freq_dict_print(countmap(insStrCol), thresh=1.1)
        end
    end

    endStrCol = [degap((alignments[j][2])[maps[j][end]+1:end]) for j in 1:length(maps)]
    if freq(endStrCol, "") < thresh
        println("After ref ends:")
        freq_dict_print(countmap(endStrCol), thresh=1.1)
    end
end

"""
    diff_in_homopolymer_region(alignment::Array{String, 1}; polylen=3)

Returns true if two aligned sequences differ only by single gaps in homopolyer regions
(ie one gap per region). `alignment` is an array of two strings that have already been aligned.
A homopolymer region is determined to be a region of a single repeated
nucleotide of length at least `polylen`.
"""
function diff_in_homopolymer_region(alignment::Array{String, 1}; polylen=3)
    for i in 1:length(alignment[1])
        if alignment[1][i] == '-' || alignment[2][i] == '-'
            lo = max(1, i - polylen + 1)
            hi = min(length(alignment[1]), i + polylen - 1)
            if !differ_by_just_one_gap([alignment[1][lo:hi], alignment[2][lo:hi]], polylen=polylen)
                return false
            end
        end
    end
    return true
end

"""
    merge_consensuses(seqs)

Takes in an array of 2-array [sequence, size] pairs and returns a list with
reach sequence unique with summed sizes among all previous copies.
"""
function merge_consensuses(seqs)
    tmp = sort(seqs, by=s->s[1]);
    i = 1
    while i < length(tmp)
        if tmp[i][1] == tmp[i+1][1]
            # consolidate cluster sizes
            tmp[i][2] += tmp[i+1][2]
            deleteat!(tmp, i+1)
        else
            i += 1
        end
    end
    return tmp
end

"""
    merge_diff_homopolymers(seqs; thresh = 0.002, k = 6, polylen = 3)

Takes in an array of 2-array [sequence, size] pairs and combines sequences that are identical
or differ by single gaps in homopolymer regions, summing respective frequencies. Slower than
naive `merge_consensuses(seqs)`.
A homopolymer region is determined to be a region of a single repeated
nucleotide of length at least `polylen`.
`thresh` gives an upper bound on distance in `k`-mer space for sequences to be compared for
homopolymer differences.
"""
function merge_diff_homopolymers(seqs; thresh = 0.002, k = 6, polylen = 3)
    # sort by decreasing size
    tmp = reverse(sort(seqs, by=s->s[2]))
    kmer_vecs = [kmer_count(s[1], k) for s in tmp]
    i = 1
    while i < length(tmp)
        j = i + 1
        while j < length(tmp)
            # if identical
            if tmp[i][1] == tmp[j][1]
                # consolidate cluster sizes
                tmp[i][2] += tmp[j][2]
                deleteat!(tmp, j)
                deleteat!(kmer_vecs, j)
            # if differ in homopolyer region
            elseif corrected_kmer_dist(kmer_vecs[i], kmer_vecs[j], k=k) < thresh &&
                    diff_in_homopolymer_region(kmer_seeded_align(tmp[i][1], tmp[j][1]), polylen=polylen)
                # consolidate cluster sizes
                tmp[i][2] += tmp[j][2]
                deleteat!(tmp, j)
                deleteat!(kmer_vecs, j)
            else
                j += 1
            end
        end
        i += 1
    end
    return tmp
end

"""
    get_coarse_centroid(seqs::Array{String, 1}; subsample = 1000, k = 4)

Returns a 'master consensus' representing largest coarse cluster of sequences,
computed among a subsample of `reads`. The consensus is the closest sequence of
`reads` to the mean `k`-mer vector in the largest computed cluster.
"""
function get_coarse_centroid(seqs::Array{String, 1}; subsample = 1000, k = 4)
    if subsample >= length(seqs)
        sampled_seqs = seqs
    else
        inds = shuffle([1:length(seqs)...])[1:subsample]
        sampled_seqs = seqs[inds]
    end
    kmer_vecs = [kmer_count(seq, k) for seq in sampled_seqs]
    @time _, sizes, _, centroid_inds = dp_centers(kmer_vecs, 0.15)
    return sampled_seqs[centroid_inds[indmax(sizes)]]
end
