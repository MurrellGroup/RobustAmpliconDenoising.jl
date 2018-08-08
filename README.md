## Synopsis

Robust amplicon denoiser (RAD) contains cluster and consensus tools. Packaging in progress, so installation is not completely straightforward yet.

## DOCS
https://murrellgroup.github.io/RAD.jl/

## Installation
We need to install 3 packages to run RAD.
```julia
Pkg.clone("https://github.com/MurrellGroup/NextGenSeqUtils.jl.git")
Pkg.clone("https://github.com/MurrellGroup/DPMeansClustering.jl.git")
Pkg.clone("https://github.com/MurrellGroup/RAD.jl.git")
```

## Load
```julia
using RAD
```

## Getting started quickly:
```julia
seqs, QVs, seq_names = read_fastq("someInputFile.fastq")
templates,template_sizes,template_indices = denoise(seqs)
write_fasta("someOutputFile.fasta",templates,names = ["seq$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])
```

## But...
You likely want to filter your reads by length, and by expected error rate. Also, PacBio reads come in random orientations, and you probably want to figure out how they should be oriented.
