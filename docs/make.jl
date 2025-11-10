using Documenter
using AnalyticEMModes

# Configurar DocMeta para el módulo
DocMeta.setdocmeta!(AnalyticEMModes, :DocTestSetup, :(using AnalyticEMModes); recursive=true)

makedocs(modules = [AnalyticEMModes],
    format = Documenter.HTML(; size_threshold=100_000_000),
    sitename = "AnalyticEMModes.jl",
    pages = Any[
        "Home" => "Introduction.md",
        "Usage" => "Usage.md",
        "Examples" => Any[
            "Gallery" => "Examples.md",
            "Rectangular" => "Rectangular/rectangular.md",
            "Cylindrical" => Any[
                "Circular Waveguides" => "Cylindrical/Circular.md",
                "Coaxial Waveguides" => "Cylindrical/Coaxial.md",
                "Radial & Wedge" => "Cylindrical/Radial.md"
            ],
            "Elliptic" => "Elliptic/elliptic.md"
        ]
    ],
    doctest = false,
    checkdocs = :none
)

deploydocs(
    repo = "github.com/uvegege/AnalyticEMModes.jl.git",
    devbranch = "main"
)
