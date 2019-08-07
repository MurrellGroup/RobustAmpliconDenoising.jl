# Robust Amplicon Denoising
This is the repository for the Robust Amplicon Denoising tool, which is suitable for denoising very long PacBio amplicon reads with relatively high error rates. Note that FAD is now available from within the NextGenSeqUtils.jl package.

## To use the RAD/FAD webserver tool, please visit:
https://biotools.ki.murrell.group/study/denoise/

## Note: This package is not yet registered for Julia v1.0. To install this package on Julia 1.0, please use the following:
```julia 
using Pkg
Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="1.0", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))
Pkg.add(PackageSpec(name="DPMeansClustering", rev="1.0", url = "https://github.com/MurrellGroup/DPMeansClustering.jl.git"))
Pkg.add(PackageSpec(name="RobustAmpliconDenoising", rev="1.0", url = "https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git"))
```

## Synopsis

Robust amplicon denoiser (RAD) contains cluster and consensus tools. Packaging in progress, so installation is not completely straightforward yet.

## DOCS
https://murrellgroup.github.io/RobustAmpliconDenoising.jl/

## Load
```julia
using RobustAmpliconDenoising
```

## Getting started quickly:
```julia
seqs, QVs, seq_names = read_fastq("someInputFile.fastq")
templates,template_sizes,template_indices = denoise(seqs)
write_fasta("someOutputFile.fasta",templates,names = ["seq$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])
```

## But...
You likely want to filter your reads by length, and by expected error rate. Also, PacBio reads come in random orientations, and you probably want to figure out how they should be oriented.

For an example for how to orient using primers, see: https://nextjournal.com/Murrell-Lab/scfv-fad-analysis/
