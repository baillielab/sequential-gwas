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
        "Prereuisistes" => [
            "ODAP" => "odap_prerequisites.md",
            "UKB RAP" => "rap_prerequisites.md",
        ],
        "Workflows" => [
            "Combining Datasets" => "combining_datasets.md",
            "Genotypes Imputation" => "genotypes_imputation.md",
            "Merging the GenOMICC and UK Biobank Cohorts" => "ukb_merge.md",
            "GWAS" => "gwas.md",
        ],
        "For Developpers" => [
            "Development" => "development.md",
            "Mock Data" => "mock_data.md",
            "Index Of Julia Functions" => "julia_fns.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/baillielab/genomicc-workflows",
    devbranch="main",
)
