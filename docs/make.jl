using Documenter
using AnalyticEMModes

# Configurar DocMeta para el módulo
DocMeta.setdocmeta!(AnalyticEMModes, :DocTestSetup, :(using AnalyticEMModes); recursive=true)

makedocs(
    clean = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://uvegege.github.io/AnalyticEMModes.jl",
        assets = String[],
        size_threshold = 100_000_000
    ),
    sitename = "AnalyticEMModes.jl",
    authors = "",
    pages = Any[
        "Home" => "Introduction.md",
        "Overview" => "Overview.md",
        "Examples" => Any[
            "Gallery" => "Ejemplos.md",
            "Cylindrical" => Any[
                "Circular Waveguides" => "Cylindrical/Circular.md",
                "Coaxial Waveguides" => "Cylindrical/Coaxial.md",
                "Radial & Wedge" => "Cylindrical/Radial.md"
            ],
            "Rectangular" => "Rectangular/rectangular.md",
            "Elliptic" => "Elliptic/elliptic.md"
        ]
    ],
    doctest = false,
    checkdocs = :none,
    remotes = nothing
)

deploydocs(
    repo = "github.com/uvegege/AnalyticEMModes.jl.git",
    devbranch = "main"
)
