using Documenter
using HydroElasticFEM

makedocs(
    sitename = "HydroElasticFEM.jl",
    authors  = "Shagun Agarwal, Oriol Colomes",
    modules  = [HydroElasticFEM],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://CMOE-TUDelft.github.io/HydroElasticFEM.jl",
        repolink   = "https://github.com/CMOE-TUDelft/HydroElasticFEM.jl",
    ),
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Guides" => [
            "Overview"                 => "guide/index.md",
            "Getting Started"          => "guide/getting_started.md",
            "First Simulation"         => "guide/first_simulation.md",
            "Examples"                 => "guide/examples.md",
            "Theory"                   => "guide/theory.md",
            "Architecture"             => "guide/architecture.md",
            "Adding a New Structure"   => "guide/adding_structure.md",
            "Debugging"                => "guide/debugging.md",
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
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/CMOE-TUDelft/HydroElasticFEM.jl.git",
        devbranch = "main",
    )
end
