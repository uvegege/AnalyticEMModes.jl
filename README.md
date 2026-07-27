# AnalyticEMModes.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://uvegege.github.io/AnalyticEMModes.jl/)


**AnalyticEMModes.jl** is a Julia package for computing analytical electromagnetic field solutions for canonical boundary value problems in multiple coordinate systems (waveguides, spherical modes, and spheroidal vector waves). All bounded geometries assume perfect electric conductor (PEC) boundary conditions.


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
| Triangular     | `tri` | (x, y) reference triangles | ✓            |
| Circular       | `cwg` | (r, θ)           | ✓                   |
| Coaxial        | `coax`| (r, θ)           | ✓                   |
| Radial         | `radial`| (ρ, φ, z)      | ✓                   |
| Wedge          | `wedge` | (ρ, φ, z)      | ✓                   |
| Elliptical       | `ewg` | (ξ, η)           | ✓                   |
| Elliptic radial | `elliptic_radial` | (ξ, η, z) elliptic cylindrical | ✓ |
| Confocal elliptic coaxial | `elliptic_coax` | (ξ, η) elliptic cylindrical annulus | ✓ |
| Spherical      | `sph` | (r, θ, φ)      | ✓                   |
| Spheroidal      | `spheroidal` | (ξ, η, ϕ) prolate/oblate | ✗ |

## Examples

Visualization examples for all supported geometries are documented in:
- [`docs/src/Examples.md`](docs/src/Examples.md)
