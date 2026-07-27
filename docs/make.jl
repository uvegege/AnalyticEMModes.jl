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
            "Triangular" => "Triangular/triangular.md",
            "Cylindrical" => Any[
                "Circular Waveguides" => "Cylindrical/Circular.md",
                "Coaxial Waveguides" => "Cylindrical/Coaxial.md",
                "Radial & Wedge" => "Cylindrical/Radial.md",
            ],
            "Elliptic" => Any[
                "Elliptic Waveguides" => "Elliptic/elliptic.md",
                "Elliptic Radial Gallery" => "Elliptic/elliptic_radial_gallery.md",
                "Confocal Elliptic Coax Gallery" => "Elliptic/confocal_elliptic_coax_gallery.md",
            ],
            "Spherical" => "Spherical/spherical.md",
            "Spheroidal" => "Spheroidal/spheroidal.md"
        ],
        "Theory" => Any[
            "Introduction" => "theory/introduction.md",
            "Coordinate Systems" => "theory/coordinate.md",
            "Mode Theory" => "theory/em_solutions.md",
            "Power Normalization" => "theory/power_normalization.md"
            ]
    ],
    doctest = false,
    checkdocs = :none
)

deploydocs(
    repo = "github.com/uvegege/AnalyticEMModes.jl.git",
    devbranch = "main",
    forcepush = true
)
