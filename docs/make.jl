using GenomiccWorkflows
using Documenter

DocMeta.setdocmeta!(GenomiccWorkflows, :DocTestSetup, :(using GenomiccWorkflows); recursive=true)

makedocs(;
    modules=[GenomiccWorkflows],
    authors="Olivier Labayle <olabayle@gmail.com> and contributors",
    sitename="GenOMICC Workflows",
    format=Documenter.HTML(;
        canonical="https://baillielab.github.io/genomicc-workflows/stable/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Walk Through" => "walk_through.md",
        "Prerequisites" => [
            "ODAP" => "odap_prerequisites.md",
            "UKB RAP" => "rap_prerequisites.md",
        ],
        "Workflows" => [
            "Combining GenOMICC Datasets" => "combining_genomicc_datasets.md",
            "GenOMICC Genotypes Imputation" => "genotypes_imputation.md",
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
