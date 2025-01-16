using SequantialGWAS
using Documenter

DocMeta.setdocmeta!(SequantialGWAS, :DocTestSetup, :(using SequantialGWAS); recursive=true)

makedocs(;
    modules=[SequantialGWAS],
    authors="Olivier Labayle <olabayle@gmail.com> and contributors",
    sitename="SequantialGWAS.jl",
    format=Documenter.HTML(;
        canonical="https://baillielab.github.io/SequantialGWAS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/baillielab/sequential-gwas",
    devbranch="main",
)
