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
        "Workflow" => "combining_genomicc_datasets.md",
        "For Developpers" => [
            "Development" => "development.md",
            "Mock Data" => "mock_data.md",
            "Legacy" => "misc.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/baillielab/genomicc-workflows",
    devbranch="main",
)
