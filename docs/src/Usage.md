# Usage Guide

This guide provides detailed information on using the AnalyticEMModes.jl package for computing electromagnetic field solutions in various geometries.

## General API Structure

For each geometry, the package provides a consistent set of functions following a standardized naming convention. The general pattern replaces `***` with the specific waveguide type abbreviation:

### Field Computation Functions

```julia
# TE mode fields (Transverse Electric: Ez = 0)
te_***_fields(coords..., geometry_params..., mode_indices..., frequency, μᵣ, εᵣ)

# TM mode fields (Transverse Magnetic: Hz = 0)
tm_***_fields(coords..., geometry_params..., mode_indices..., frequency, μᵣ, εᵣ)

```

### Cutoff Wavenumber Functions

```julia
# Cutoff wavenumber for TE/TM modes
kc = kc_***(geometry_params..., mode_indices..., mode_type)
```

### Mode Sorting Functions

```julia
# Get first N modes ordered by cutoff frequency
modes = first_n_modes_***(N, geometry_params...)
```

## Coordinate Systems and Parameters

### Rectangular Waveguides (`rwg`)

**Coordinates:** Cartesian (x, y)
- `x`: Position along width (0 ≤ x ≤ W)
- `y`: Position along height (0 ≤ y ≤ H)

**Geometry Parameters:**
- `W`: Waveguide width (m)
- `H`: Waveguide height (m)

**Mode Indices:**
- `m`: Number of half-wave variations along width
- `n`: Number of half-wave variations along height

**Functions:**
```julia
Ex, Ey, Ez, Hx, Hy, Hz = te_rwg_fields(x, y, W, H, m, n, f, μᵣ, εᵣ)
Ex, Ey, Ez, Hx, Hy, Hz = tm_rwg_fields(x, y, W, H, m, n, f, μᵣ, εᵣ)
kc = kc_rwg(W, H, m, n)
modes = first_n_modes_rwg(N, W, H)
```

### Circular Waveguides (`cwg`)

**Coordinates:** Cylindrical (r, θ)
- `r`: Radial distance from center (0 ≤ r ≤ R)
- `θ`: Azimuthal angle (radians)

**Geometry Parameters:**
- `R`: Waveguide radius (m)

**Mode Indices:**
- `m`: Azimuthal mode number (angular variations)
- `n`: Radial mode number (radial zeros)

**Functions:**
```julia
Er, Eθ, Ez, Hr, Hθ, Hz = te_cwg_fields(r, θ, R, m, n, f, μᵣ, εᵣ)
Er, Eθ, Ez, Hr, Hθ, Hz = tm_cwg_fields(r, θ, R, m, n, f, μᵣ, εᵣ)
kc = kc_cwg(R, m, n, mode_type)  # mode_type = :TE or :TM
modes = first_n_modes_cwg(N, R)
```

### Coaxial Waveguides (`coax`)

**Coordinates:** Cylindrical (r, θ)
- `r`: Radial distance from center (R2 ≤ r ≤ R1)
- `θ`: Azimuthal angle (radians)

**Geometry Parameters:**
- `R1`: Outer conductor radius (m)
- `R2`: Inner conductor radius (m)

**Mode Indices:**
- `m`: Azimuthal mode number
- `n`: Radial mode number

**Functions:**
```julia
Er, Eθ, Ez, Hr, Hθ, Hz = te_coax_fields(r, θ, R1, R2, m, n, f, μᵣ, εᵣ)
Er, Eθ, Ez, Hr, Hθ, Hz = tm_coax_fields(r, θ, R1, R2, m, n, f, μᵣ, εᵣ)
Er, Eθ, Ez, Hr, Hθ, Hz = tem_coax_fields(r, θ, R1, R2, f, μᵣ, εᵣ)
kc = kc_coax(R1, R2, m, n, mode_type) # mode_type = :TE or :TM
modes = first_n_modes_coax(N, R1, R2)
```

**Special:** TEM mode has no cutoff (kc = 0).

### Elliptical Waveguides (`ewg`)

**Coordinates:** Elliptic cylindrical (ξ, η)
- `ξ`: Radial-like elliptic coordinate
- `η`: Angular-like elliptic coordinate

**Geometry Parameters:**
- `a`: Semi-major axis (m)
- `b`: Semi-minor axis (m)

**Mode Indices:**
- `m`: Angular Mathieu function order
- `n`: Radial Mathieu function order
- `even`: Boolean (true = even modes, false = odd modes)

**Functions:**
```julia
Eξ, Eη, Ez, Hξ, Hη, Hz = te_ewg_fields(ξ, η, a, b, m, n, even, f, μᵣ, εᵣ)
Eξ, Eη, Ez, Hξ, Hη, Hz = tm_ewg_fields(ξ, η, a, b, m, n, even, f, μᵣ, εᵣ)
kc = kc_ewg(a, b, m, n, even, mode_type)
modes = first_n_modes_ewg(N, a, b)
```

**Coordinate Conversion:**
```julia
ξ, η = cart2elliptic(x, y, a, b)
```

### Radial Waveguides (`radial`)

**Coordinates:** Cylindrical (ρ, φ, z)
- `ρ`: Radial distance (propagation direction)
- `φ`: Azimuthal angle
- `z`: Vertical position (0 ≤ z ≤ H)

**Geometry Parameters:**
- `H`: Height between parallel plates (m)
- `m`: Azimuthal mode number
- `n`: Vertical mode number
- `ρ0`: Reference radial position

**Functions:**
```julia
Eρ, Eφ, Ez, Hρ, Hφ, Hz = te_radial_fields(ρ, φ, z, H, m, n, ρ0, kρ, f, μᵣ, εᵣ)
Eρ, Eφ, Ez, Hρ, Hφ, Hz = tm_radial_fields(ρ, φ, z, H, m, n, ρ0, kρ, f, μᵣ, εᵣ)
```

### Wedge Waveguides (`wedge`)

**Coordinates:** Cylindrical sector (ρ, φ, z)
- Angular extent: 0 ≤ φ ≤ φ₀

**Geometry Parameters:**
- `H`: Height between parallel plates (m)
- `φ0`: Wedge angle (radians)
- `p`: Angular mode number (quantized by φ₀)
- `n`: Vertical mode number

**Functions:**
```julia
Eρ, Eφ, Ez, Hρ, Hφ, Hz = te_wedge_fields(ρ, φ, z, H, φ0, p, n, kρ, ρ0, f, μᵣ, εᵣ)
Eρ, Eφ, Ez, Hρ, Hφ, Hz = tm_wedge_fields(ρ, φ, z, H, φ0, p, n, kρ, ρ0, f, μᵣ, εᵣ)
```


## Return Values

Field functions return a 6-tuple of complex-valued phasor components:

```julia
(E_coord1, E_coord2, E_coord3, H_coord1, H_coord2, H_coord3)
```

The coordinate system depends on the waveguide geometry:
- **Cartesian:** (Ex, Ey, Ez, Hx, Hy, Hz)
- **Cylindrical:** (Er, Eθ, Ez, Hr, Hθ, Hz)
- **Elliptic:** (Eξ, Eη, Ez, Hξ, Hη, Hz)

## Helper Functions

### Cutoff Frequency and Propagation

```julia
using AnalyticEMModes

# Convert cutoff wavenumber to cutoff frequency
fc = cutoff_frequency(kc, μᵣ, εᵣ)

# Compute propagation constant
β = phase_constant(kc, f, μᵣ, εᵣ)

```

### Coordinate Transformations

```julia
# Cylindrical unit vectors and metric
h_r, h_θ, e_r_x, e_r_y, e_θ_x, e_θ_y = metric_and_unit_cylindrical(r, θ)

# Elliptic unit vectors and metric
h, e_ξ_x, e_ξ_y, e_η_x, e_η_y = metric_and_unit_elliptic(ρ, ξ, η)

# Convert field components from cylindrical to Cartesian
Fx = Fr * e_r_x + Fθ * e_θ_x
Fy = Fr * e_r_y + Fθ * e_θ_y
```

## Mode Sorting and Selection

Finding modes ordered by cutoff frequency:

```julia
# Rectangular waveguide
W, H = 0.05, 0.025
modes = first_n_modes_rwg(10, W, H)
for (mode_type, m, n, kc) in modes
    fc = cutoff_frequency(kc, 1.0, 1.0)
    println("$mode_type($m,$n): fc = $(fc/1e9) GHz")
end
```

### Evaluate Fields

```julia
# Single point
r, θ = 0.005, π/4
Er, Eθ, Ez, Hr, Hθ, Hz = te_cwg_fields(r, θ, R, m, n, f, μᵣ, εᵣ)

# Multiple points (vectorized)
r_vec = range(0, R, length=100)
θ_vec = fill(π/4, 100)
fields = te_cwg_fields(r_vec, θ_vec, R, m, n, f, μᵣ, εᵣ)
```

## Examples

For detailed examples with visualizations, see:
- [Circular Waveguides](Cylindrical/Circular.md)
- [Coaxial Waveguides](Cylindrical/Coaxial.md)
- [Radial & Wedge](Cylindrical/Radial.md)
- [Rectangular Waveguides](Rectangular/rectangular.md)
- [Elliptic Waveguides](Elliptic/elliptic.md)


