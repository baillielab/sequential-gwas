using SequentialGWAS
using Documenter

DocMeta.setdocmeta!(SequentialGWAS, :DocTestSetup, :(using SequentialGWAS); recursive=true)

makedocs(;
    modules=[SequentialGWAS],
    authors="Olivier Labayle <olabayle@gmail.com> and contributors",
    sitename="SequentialGWAS.jl",
    format=Documenter.HTML(;
        canonical="https://baillielab.github.io/SequentialGWAS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Cohorts Description" => "cohorts.md",
        "Input Data" => "input_data.md",
        "Run" => "run.md",
        "Outputs" => "outputs.md",
        "Development" => "development.md",
        "Mock Data" => "mock_data.md",
    ],
)

deploydocs(;
    repo="github.com/baillielab/sequential-gwas",
    devbranch="main",
)
