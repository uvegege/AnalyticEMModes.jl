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
| Circular       | `cwg` | (r, θ)           | ✓                   |
| Coaxial        | `coax`| (r, θ)           | ✓                   |
| Radial         | `radial`| (ρ, φ, z)      | ✓                   |
| Wedge          | `wedge` | (ρ, φ, z)      | ✓                   |
| Elliptical       | `ewg` | (ξ, η)           | ✓                   |
| Spherical      | `sph` | (r, θ, φ)      | ✓                   |
| Spheroidal      | `spheroidal` | (ξ, η, ϕ) prolate/oblate | ✗ |

## Spheroidal API (current scope)

The spheroidal module currently exposes vector-wave building blocks:

- `mn_spheroidal_vector`, `m_spheroidal_vector`, `n_spheroidal_vector`
- `mn_spheroidal_vectors`, `m_spheroidal_vectors`, `n_spheroidal_vectors`
- `mn_spheroidal_vectors_mnmax`
- `ProlateSpheroidalBasis`, `OblateSpheroidalBasis`, `SpheroidalB`
- coordinate helpers (`pro2cart`, `obl2cart`, `cart2pro`, `cart2obl`, scale factors)


## Examples

Visualization examples for all supported geometries are documented in:
- [`docs/src/Examples.md`](docs/src/Examples.md)
