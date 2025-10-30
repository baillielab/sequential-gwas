using Documenter

makedocs(;
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
        "ODAP" => "odap_prerequisites.md",
        "Workflows" => [
            "Combining GenOMICC Datasets" => "combining_genomicc_datasets.md",
            "GenOMICC Genotypes Imputation" => "genotypes_imputation.md",
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
