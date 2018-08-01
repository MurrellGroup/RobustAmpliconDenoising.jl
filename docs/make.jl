# Run `julia make.jl` in this folder to generate .html pages in a build/
# directory, then open index.md for docs

push!(LOAD_PATH,"../src/")
using Documenter, RAD

makedocs(
    format = :html,
    sitename = "RAD.jl",
    modules = [RAD],
    pages = [
        "index.md",
        "Documentation" => [
            "clusterpipeline.md",
            "consensus.md"
        ]
    ]
)
