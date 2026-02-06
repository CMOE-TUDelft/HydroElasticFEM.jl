using HydroElasticFEM
using Documenter

makedocs(;
    modules=[HydroElasticFEM],
    authors="Shagun Agarwal <shagun.1994@gmail.com>, Oriol Colomes <J.O.ColomesGene@tudelft.nl> and contributors",
    repo="https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/blob/{commit}{path}#{line}",
    sitename="HydroElasticFEM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://CMOE-TUDelft.github.io/HydroElasticFEM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/CMOE-TUDelft/HydroElasticFEM.jl",
    devbranch="main",
)
