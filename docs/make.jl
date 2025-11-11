using Documenter
using AnalyticEMModes

# Configurar DocMeta para el módulo
DocMeta.setdocmeta!(AnalyticEMModes, :DocTestSetup, :(using AnalyticEMModes); recursive=true)

makedocs(modules = [AnalyticEMModes],
    format = Documenter.HTML(; size_threshold=100_000_000),
    clean = true,
    sitename = "AnalyticEMModes.jl",
    pages = Any[
        "index.md",
        "Usage.md",
        "api_reference.md",
        "Examples" => Any[
            "Gallery" => "Examples.md",
            "Rectangular" => "Rectangular/rectangular.md",
            "Cylindrical" => Any[
                "Circular Waveguides" => "Cylindrical/Circular.md",
                "Coaxial Waveguides" => "Cylindrical/Coaxial.md",
                "Radial & Wedge" => "Cylindrical/Radial.md"
            ],
            "Elliptic" => "Elliptic/elliptic.md"
        ],
        "Theory" => Any[
            "Introduction" => "theory/introduction.md",
            "Coordinate Systems" => "theory/coordinate.md",
            "Mode Theory" => "theory/em_solutions.md",
            "Rectangular" => "theory/theory_rectangular.md"
            ]
    ],
    doctest = false,
    checkdocs = :none
)

deploydocs(
    repo = "github.com/uvegege/AnalyticEMModes.jl.git",
    devbranch = "main"
)
