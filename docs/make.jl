using Documenter, RAD
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/USER_NAME/PACKAGE_NAME.jl.git",
    julia  = "nightly",
    osname = "osx"
makedocs()
