using Documenter
using HydroElasticFEM

makedocs(
    sitename = "HydroElasticFEM.jl",
    authors  = "Shagun Agarwal, Oriol Colomes",
    modules  = [HydroElasticFEM],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://github.com/CMOE/HydroElasticFEM.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Guides" => [
            "Getting Started" => "guide/getting_started.md",
            "Examples"        => "guide/examples.md",
            "Theory"          => "guide/theory.md",
        ],
        "API Reference" => [
            "Overview"          => "api/index.md",
            "Physics"           => "api/physics.md",
            "Geometry"          => "api/geometry.md",
            "Simulation"        => "api/simulation.md",
            "Assembly Contexts" => "api/assembly_contexts.md",
            "Parameter Handler" => "api/parameter_handler.md",
        ],
    ],
    checkdocs = :exports,
    doctest   = true,
    strict    = false,
)

deploydocs(
    repo = "github.com/CMOE/HydroElasticFEM.jl.git",
    devbranch = "main",
)
