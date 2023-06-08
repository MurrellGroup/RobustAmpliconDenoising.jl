# Robust Amplicon Denoising
This is the repository for the Robust Amplicon Denoising tool, which is suitable for denoising very long PacBio amplicon reads with relatively high error rates. Note that FAD is now available from within the NextGenSeqUtils.jl package.

## To use the RAD/FAD webserver tool, please visit:
https://biotools.ki.murrell.group/study/denoise/

## Note: This package is not registered. To install this package, along with its dependencies, please use the following:
```julia 
using Pkg
Pkg.add(url="https://github.com/MurrellGroup/NextGenSeqUtils.jl")
Pkg.add(url="https://github.com/MurrellGroup/DPMeansClustering.jl")
Pkg.add(url="https://github.com/MurrellGroup/RobustAmpliconDenoising.jl")
```

## Getting started quickly:
Load the package:
```julia
using NextGenSeqUtils, RobustAmpliconDenoising
```
The run, via:
```julia
seqs, QVs, seq_names = read_fastq("someInputFile.fastq")
templates,template_sizes,template_indices = denoise(seqs)
write_fasta("someOutputFile.fasta",templates,names = ["seq$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])
```

## But...
You likely want to filter your reads by length, and by expected error rate. Also, PacBio reads come in random orientations, and you probably want to figure out how they should be oriented.

For an example for how to orient using primers, see: https://nextjournal.com/Murrell-Lab/scfv-fad-analysis/

## DOCS
https://murrellgroup.github.io/RobustAmpliconDenoising.jl/
