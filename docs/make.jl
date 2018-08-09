using Documenter, RAD
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/MurrellGroup/RAD.jl.git",
    julia  = "nightly",
    osname = "osx"
makedocs()
