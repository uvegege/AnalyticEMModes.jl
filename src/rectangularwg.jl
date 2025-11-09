"""
    kc_rwg(a, b, m, n)

Calculates the cutoff wave number of the m,n mode for a rectangular waveguide with dimensions a, b where a >= b.

# Arguments
- `a`: Width of the waveguide (larger dimension)
- `b`: Height of the waveguide (smaller dimension)
- `m`: Mode index in x-direction
- `n`: Mode index in y-direction
"""
kc_rwg(a, b, m, n) = π * hypot(m/a, n/b)


"""
    te_rwg_modal_f(x, y, a, b, m, n)

Modal functions of a TE_{m,n} mode for a rectangular guide with dimensions a, b with a >= b.

# Arguments
- `x`, `y`: Coordinates where x ∈ [0, a] and y ∈ [0, b]
- `a`: Width of the waveguide
- `b`: Height of the waveguide
- `m`, `n`: Mode indices

# Returns
A tuple (∂Fz/∂dx, ∂Fz/∂dy, Fz) for a given (x, y).
"""
function te_rwg_modal_f(x, y, a, b, m::Integer, n::Integer)
    ∂ψᵢ = -sinpi(m*x/a) * cospi(n*y/b) * m*π/a # ∂Fz/∂dx
    ∂ψⱼ = -cospi(m*x/a) * sinpi(n*y/b) * n*π/b # ∂Fz/∂dy
    ψₖ =  cospi(m*x/a) * cospi(n*y/b) # Fz
    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end

"""
    tm_rwg_modal_f(x, y, a, b, m, n)

Modal functions of a TM_{m,n} mode for a rectangular guide with dimensions a, b with a >= b.

# Arguments
- `x`, `y`: Coordinates where x ∈ [0, a] and y ∈ [0, b]
- `a`: Width of the waveguide
- `b`: Height of the waveguide
- `m`, `n`: Mode indices

# Returns
A tuple (∂Fz/∂dx, ∂Fz/∂dy, Fz) for a given (x, y).
"""
function tm_rwg_modal_f(x, y, a, b, m::Integer, n::Integer)
    ∂ψᵢ = cospi(m*x/a) * sinpi(n*y/b) * m*π/a # ∂Fz/∂dx
    ∂ψⱼ = cospi(n*y/b) * sinpi(m*x/a) * n*π/b # ∂Fz/∂dy
    ψₖ = sinpi(m*x/a) * sinpi(n*y/b) # Fz
    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end


"""
    te_rwg_fields(x, y, a, b, m, n, c_e, c_h)
    te_rwg_fields(x, y, a, b, m, n, f, μᵣ, εᵣ)

Calculate the electric and magnetic field components of a TE_{m,n} mode.

# Arguments
- `x`, `y`: Coordinates where x ∈ [0, a] and y ∈ [0, b]
- `a`: Width of the waveguide
- `b`: Height of the waveguide
- `m`, `n`: Mode indices
- `c_e`, `c_h`: Field coupling coefficients (from `te_coefficients`), or
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Returns
`(Ex, Ey, Ez, Hx, Hy, Hz)`
"""
function te_rwg_fields(x, y, a, b, m, n, c_e, c_h)
   
    ∂ψᵢ, ∂ψⱼ, ψₖ = te_rwg_modal_f(x, y, a, b, m, n)

    Ex = -c_e * ∂ψᵢ
    Ey = +c_e * ∂ψⱼ
    Ez = zero(Ex)
    Hx = -c_h * ∂ψⱼ
    Hy = -c_h * ∂ψᵢ
    Hz =  -im *  ψₖ

    return (Ex, Ey, Ez, Hx, Hy, Hz)
end

function te_rwg_fields(x, y, a, b, m, n, f, μᵣ, εᵣ)
    
    kc = kc_rwg(a, b, m, n)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)
    (Ex, Ey, Ez, Hx, Hy, Hz) = te_rwg_fields(x, y, a, b, m, n, c_e, c_h)

    return (Ex, Ey, Ez, Hx, Hy, Hz) 
end

function te_rwg_fields(x::AbstractArray{T, N}, y::AbstractArray{T, N}, a, b, m, n, f, μᵣ, εᵣ) where {T, N}
    
    kc = kc_rwg(a, b, m, n)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(x, NTuple{6, Complex{T}})
    for idx in eachindex(x)
        fields[idx] = te_rwg_fields(x[idx], y[idx], a, b, m, n, c_e, c_h)
    end
    return fields
end

"""
    tm_rwg_fields(x, y, a, b, m, n, c_e, c_h)
    tm_rwg_fields(x, y, a, b, m, n, f, μᵣ, εᵣ)

Compute the electric and magnetic field components of a TM_{m,n} mode.

# Arguments
- `x`, `y`: Coordinates where x ∈ [0, a] and y ∈ [0, b]
- `a`: Width of the waveguide
- `b`: Height of the waveguide
- `m`, `n`: Mode indices
- `c_e`, `c_h`: Field coupling coefficients (from `tm_coefficients`), or
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity

# Returns
`(Ex, Ey, Ez, Hx, Hy, Hz)`

# Example
```julia
c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)
tm_rwg_fields(x, y, a, b, m, n, c_e, c_h)
```
"""
function tm_rwg_fields(x, y, a, b, m, n, c_e, c_h)
    
    ∂ψᵢ, ∂ψⱼ, ψₖ = tm_rwg_modal_f(x, y, a, b, m, n)
    
    Ex = -c_e * ∂ψᵢ
    Ey = -c_e * ∂ψⱼ
    Ez =  -im *  ψₖ
    Hx = +c_h * ∂ψⱼ
    Hy = -c_h * ∂ψᵢ
    Hz = zero(Ex)
    
    return (Ex, Ey, Ez, Hx, Hy, Hz)
end

function tm_rwg_fields(x, y, a, b, m, n, f, μᵣ, εᵣ)
    
    kc = kc_rwg(a, b, m, n)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)
    (Ex, Ey, Ez, Hx, Hy, Hz) = tm_rwg_fields(x, y, a, b, m, n, c_e, c_h)

    return (Ex, Ey, Ez, Hx, Hy, Hz) 
end

function tm_rwg_fields(x::AbstractArray{T, N}, y::AbstractArray{T, N}, a, b, m, n, f, μᵣ, εᵣ) where {T, N}
    
    kc = kc_rwg(a, b, m, n)
    β = phase_constant(kc, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(kc, β, f, μᵣ, εᵣ)

    fields = similar(x, NTuple{6, Complex{T}})
    for idx in eachindex(x)
        fields[idx] = tm_rwg_fields(x[idx], y[idx], a, b, m, n, c_e, c_h)
    end
    return fields
end

"""
    te_normalization_rwg(a, b, m, n, kc, β, f, μᵣ, εᵣ)

Normalization factor for TE modes to achieve unit power.

The normalization is derived from integrating the Poynting vector over
the waveguide cross-section, giving the total transmitted power.
In TE modes, Hₒ is normalized. δm, δn account for mode indices being 0 or non-zero.

# Arguments
- `a`: Width of the waveguide
- `b`: Height of the waveguide
- `m`, `n`: Mode indices
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function te_normalization_rwg(a, b, m, n, kc, β, f, μᵣ, εᵣ)
    ω = 2π*f
    μ = μᵣ * _μₒ
    δm = ifelse(n >= 1, 1, 2)
    δn = ifelse(m >= 1, 1, 2)
    Hₒ = sqrt( (8*kc^4) / (β*ω*μ*π^2 * (δm * m^2*b/a + δn*n^2*a/b)) )
    return Hₒ
end


"""
    tm_normalization_rwg(a, b, m, n, kc, β, f, μᵣ, εᵣ)

Normalization factor for TM modes to achieve unit power.

The normalization is derived from integrating the Poynting vector over
the waveguide cross-section, giving the total transmitted power.
In TM modes, Eₒ is normalized. In this case there is no δm, δn because m and n >= 1.

# Arguments
- `a`: Width of the waveguide
- `b`: Height of the waveguide
- `m`, `n`: Mode indices
- `kc`: Cutoff wavenumber
- `β`: Phase constant
- `f`: Frequency in Hz
- `μᵣ`: Relative permeability
- `εᵣ`: Relative permittivity
"""
function tm_normalization_rwg(a, b, m, n, kc, β, f, μᵣ, εᵣ)
    ω = 2π*f
    ε = εᵣ * _εₒ
    μ = μᵣ * _μₒ
    Eₒ = sqrt( (8*kc^4) / (β*ω*ε*π^2 * (m^2*b/a + n^2*a/b)) )
    return Eₒ
end