# AnalyticEMModes.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://uvegege.github.io/AnalyticEMModes.jl/)


**AnalyticEMModes.jl** is a Julia package for computing analytical electromagnetic field solutions for canonical boundary value problems in various geometries (such as waveguides, spherical modes, or radial propagation). All bounded geometries assume perfect electric conductor (PEC) boundary conditions.


## Installation

From the Julia REPL:

```julia
using Pkg
Pkg.add("AnalyticEMModes")
```

Or in package mode (type `]`):

```
pkg> add AnalyticEMModes
```

## Supported Geometries

| Type | Alias | Coordinate System | Power Normalization |
|:---------------|:------|:------------------|:-------------------:|
| Rectangular    | `rwg` | (x, y)           | ✓                   |
| Circular       | `cwg` | (r, θ)           | ✓                   |
| Coaxial        | `coax`| (r, θ)           | ✓                   |
| Radial         | `radial`| (ρ, φ, z)      | ✓                   |
| Wedge          | `wedge` | (ρ, φ, z)      | ✓                   |
| Elliptical       | `ewg` | (ξ, η)           | ✗                   |
| Spherical      | `sph` | (r, θ, φ)      | ✓                   |
| Spheroidal      | `spheroidal` | (_, _, _)      | 🚧 Under development |