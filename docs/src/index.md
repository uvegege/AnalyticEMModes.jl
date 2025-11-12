# AnalyticEMModes.jl

**AnalyticEMModes.jl** is a Julia package for computing analytical electromagnetic field solutions for canonical boundary value problems in various geometries (such as waveguides, spherical modes, or radial propagation). All bounded geometries assume perfect electric conductor (PEC) boundary conditions.

## Installation

From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/uvegege/AnalyticEMModes.jl")
```

Or in package mode (type `]`):

```
pkg> add https://github.com/uvegege/AnalyticEMModes.jl
```

## Overview

This package provides closed-form solutions for electric and magnetic field distributions in waveguides, enabling fast and accurate computation of mode patterns, cutoff frequencies, and propagation characteristics. Unlike numerical methods (FEM, FDTD), analytical solutions offer:

- **Exact results** (within numerical precision of special function evaluation)
- **Physical insight** through explicit functional forms
- **Benchmarking** for validating numerical solvers
- **Integration with numerical methods**, e.g. using analytical modes to excite ports, impose boundary conditions, etc.

## Supported Geometries

| Waveguide Type | Alias | Coordinate System | Power Normalization |
|:---------------|:------|:------------------|:-------------------:|
| Rectangular    | `rwg` | (x, y)           | Yes                 |
| Circular       | `cwg` | (r, θ)           | Yes                 |
| Coaxial        | `coax`| (r, θ)           | Yes                 |
| Radial         | `radial`| (ρ, φ, z)      | Yes                 |
| Wedge          | `wedge` | (ρ, φ, z)      | Yes                 |
| Elliptical       | `ewg` | (ξ, η)           | No                  |
| Spherical      | —     | (r, θ, φ)        | Under development   |


## Features

- **Field Components**: Complete 6-component electromagnetic field vectors (Er, Eθ, Ez, Hr, Hθ, Hz)
- **Cutoff Analysis**: Cutoff wavenumbers, frequencies, and wavelengths
- **Propagation Parameters**: Phase constant, attenuation, impedance
- **Mode Sorting**: Efficient algorithms to find first N modes ordered by cutoff frequency
- **Coordinate Transformations**: Conversion between Cartesian and curvilinear coordinates
- **High Accuracy**: Uses `Bessels.jl` for reliable special function evaluation


## Quick Start

### Example 1: TE₁₁ Mode in Circular Waveguide

```julia
using AnalyticEMModes

# Waveguide parameters
a = 0.01         # Radius in meters
m, n = 1, 1      # Mode indices (TE₁₁)
f = 10e9         # Frequency in Hz
μᵣ, εᵣ = 1.0, 1.0  # Relative permeability and permittivity

# Compute cutoff frequency
kc = kc_cwg(a, m, n, :TE)
fc = cutoff_frequency(kc, μᵣ, εᵣ)
println("TE₁₁ cutoff frequency: $(fc/1e9) GHz")

# Evaluate fields at point (r, θ)
r, θ = 0.005, π/4
Er, Eθ, Ez, Hr, Hθ, Hz = te_cwg_fields(r, θ, a, m, n, f, μᵣ, εᵣ)

println("Magnetic field Hz: $Hz")
```

### Example 2: Finding First N Modes

```julia
using AnalyticEMModes

# Find first 10 modes in rectangular waveguide
a, b = 0.02286, 0.01016  # WR-90 dimensions (m)
modes = first_n_modes_rwg(10, a, b)
for (i, (kind, m, n, kc)) in enumerate(modes)
    fc = cutoff_frequency(kc, 1.0, 1.0)
    println("Mode $i: $kind($m,$n), fc = $(fc/1e9) GHz")
end
```

## Documentation Structure

This documentation is organized as follows:

### User Guide
- **[Usage](Usage.md)**: Detailed function reference and conventions
- **[Examples](Examples.md)**: Gallery of mode visualizations and comparisons

### Examples by Geometry
- **[Circular Waveguides](Cylindrical/Circular.md)**: TE and TM modes in circular guides
- **[Coaxial Waveguides](Cylindrical/Coaxial.md)**: Including TEM mode
- **[Radial & Wedge](Cylindrical/Radial.md)**: Radial propagation geometries
- **[Rectangular Waveguides](Rectangular/rectangular.md)**: Standard rectangular guides
- **[Elliptic Waveguides](Elliptic/elliptic.md)**: Mathieu function solutions

### Theory (Recommended Reading)
For users interested in the mathematical foundations:
- Coordinate systems and transformations
- Power normalization 
- [WIP] Potential formulations (Hertzian potentials)  , Boundary value problems and separation of variables  , Special functions (Bessel, Mathieu) , References to classic textbooks 

## API Conventions

### Naming Conventions
- Functions use lowercase with underscores: `te_cwg_fields`, `kc_rwg`
- Mode types are symbols: `:TE`, `:TM`, `:TEM`
- Coordinate parameters follow physics notation: r, θ, z for cylindrical

### Return Values
- Field functions return 6-tuple: `(Er, Eθ, Ez, Hr, Hθ, Hz)` or `(Ex, Ey, Ez, Hx, Hy, Hz)`
- All fields are complex-valued (phasor notation)
- Units follow SI: meters, Hz, A/m, V/m

### Mode Indexing
- `(m, n)`: m = azimuthal/angular index, n = radial/vertical index
- Mode types with polarization:
  - `:TEc`, `:TEs`: Cosine/sine polarizations (circular/coaxial, same cutoff)
  - `:TEe`, `:TEo`: Even/odd modes (elliptic, different cutoffs)


## Validation

Analytical results can be validated against numerical solvers.

Example comparison code is provided in the documentation for each geometry.

## Contributing

Contributions are welcome! Areas for expansion:
- Additional geometries (spherical, etc.) with analytical solutions.
- Additional visualization utilities, like field line plotting.
- And anything else you think might fit in AnalyticEMMode.jl


## Acknowledgments

This package builds upon:
- `Bessels.jl` for Bessel function evaluation
- Classic EM and waveguide books.

## Support

- **Documentation**: This site
- **Issues**: [GitHub Issues](https://github.com/uvegege/AnalyticEMModes.jl/issues)
- **Discussions**: [GitHub Discussions](https://github.com/uvegege/AnalyticEMModes.jl/discussions)
