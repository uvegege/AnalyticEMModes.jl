## Spheroidal coordinates functions

"""
    obl2cart(a, ξ, η, ϕ)

Convert oblate spheroidal coordinates `(ξ, η, ϕ)` to Cartesian `(x, y, z)`.
`a` is the focal distance. Ranges: `ξ ≥ 0`, `η ∈ [-1, 1]`, `ϕ ∈ [0, 2π)`.
"""
function obl2cart(a, ξ, η, ϕ)
    ρ = a*sqrt(ξ^2 + 1)*sqrt(1 - η^2)
    x = ρ * cos(ϕ)
    y = ρ * sin(ϕ)
    z = a*ξ*η
    return x,y,z
end

"""
    pro2cart(a, ξ, η, ϕ)

Convert prolate spheroidal coordinates `(ξ, η, ϕ)` to Cartesian `(x, y, z)`.
`a` is the focal distance. Ranges: `ξ ≥ 1`, `η ∈ [-1, 1]`, `ϕ ∈ [0, 2π)`.
"""
function pro2cart(a, ξ, η, ϕ)
    ρ = a * sqrt(ξ^2 - 1) * sqrt(1 - η^2)
    x = ρ * cos(ϕ)
    y = ρ * sin(ϕ)  
    z = a * ξ * η
    return x, y, z
end

"""
    cart2pro(a, x, y, z)

Convert Cartesian coordinates `(x, y, z)` to prolate spheroidal `(ξ, η, ϕ)`.
`a` is the focal distance.
"""
function cart2pro(a, x, y, z)
    ρ = sqrt(x^2 + y^2)
    ϕ = atan(y, x)
    r1 = sqrt(ρ^2 + (z - a)^2)
    r2 = sqrt(ρ^2 + (z + a)^2)
    ξ = (r1 + r2) / (2*a)
    η = (r2 - r1) / (2*a)
    ξ = max(ξ, 1.0)
    η = clamp(η, -1.0, 1.0)
    return ξ, η, ϕ
end

"""
    cart2obl(a, x, y, z)

Convert Cartesian coordinates `(x, y, z)` to oblate spheroidal `(ξ, η, ϕ)`.
`a` is the focal distance.
"""
function cart2obl(a, x, y, z)
    ρ = sqrt(x^2 + y^2)
    ϕ = atan(y, x)
    
    r = sqrt(x^2 + y^2 + z^2)
    ξ = sqrt((r^2 - a^2 + sqrt((r^2 - a^2)^2 + 4*a^2*z^2)) / (2*a^2))
    η = z / (a * ξ)
     
    ξ = max(ξ, 0.0)  # ξ ≥ 0
    η = clamp(η, -1.0, 1.0)  # η ∈ [-1, 1]
    
    return ξ, η, ϕ
end

"""
    spheroidal_parameter(major_axis, minor_axis)

Compute the focal distance `a` of a spheroid from its full major and minor axes lengths.
"""
function spheroidal_parameter(major_axis, minor_axis)
    a = major_axis / 2
    b = minor_axis / 2
    focal_distance = 2 * sqrt(a^2 - b^2)  
    return focal_distance / 2  
end

"""
    scale_factors_prolate(a, ξ, η)

Return the metric scale factors `(hξ, hη, hϕ)` for prolate spheroidal coordinates at `(ξ, η)`.
"""
function scale_factors_prolate(a, ξ, η)
    hξ = a * sqrt((ξ^2 - η^2) / (ξ^2 - 1))
    hη = a * sqrt((ξ^2 - η^2) / (1 - η^2))
    hϕ = a * sqrt((ξ^2 - 1) * (1 - η^2))
    return hξ, hη, hϕ
end

"""
    scale_factors_oblate(a, ξ, η)

Return the metric scale factors `(hξ, hη, hϕ)` for oblate spheroidal coordinates at `(ξ, η)`.
"""
function scale_factors_oblate(a, ξ, η)
    hξ = a * sqrt((ξ^2 + η^2) / (ξ^2 + 1))
    hη = a * sqrt((ξ^2 + η^2) / (1 - η^2))
    hϕ = a * sqrt((ξ^2 + 1) * (1 - η^2))
    return hξ, hη, hϕ
end

function prolate_partials(a, ξ, η, ϕ)
    R = a*sqrt((ξ^2-1)*(1-η^2))
    dR_dξ = a*(ξ*(1-η^2))/sqrt((ξ^2-1)*(1-η^2))
    dR_dη = -a*(η*(ξ^2-1))/sqrt((ξ^2-1)*(1-η^2))

    ex_x =  (dR_dξ*cos(ϕ))
    ex_y =  (dR_dξ*sin(ϕ))
    ex_z =  a*η

    eη_x =  (dR_dη*cos(ϕ))
    eη_y =  (dR_dη*sin(ϕ))
    eη_z =  a*ξ

    eϕ_x = -R*sin(ϕ)
    eϕ_y =  R*cos(ϕ)
    eϕ_z =  0.0

    return ex_x, ex_y, ex_z, eη_x, eη_y, eη_z, eϕ_x, eϕ_y, eϕ_z
end

function oblate_partials(a, ξ, η, ϕ)
    R = a*sqrt((ξ^2+1)*(1-η^2))
    dR_dξ = a*(ξ*(1-η^2))/sqrt((ξ^2+1)*(1-η^2))
    dR_dη = -a*(η*(ξ^2+1))/sqrt((ξ^2+1)*(1-η^2))

    ex_x = dR_dξ*cos(ϕ)
    ex_y = dR_dξ*sin(ϕ)
    ex_z = a*η

    eη_x = dR_dη*cos(ϕ)
    eη_y = dR_dη*sin(ϕ)
    eη_z = a*ξ

    eϕ_x = -R*sin(ϕ)
    eϕ_y =  R*cos(ϕ)
    eϕ_z =  0.0

    return ex_x, ex_y, ex_z,
           eη_x, eη_y, eη_z,
           eϕ_x, eϕ_y, eϕ_z
end

"""
    prolate_vector_to_cartesian(a, ξ, η, ϕ, Vξ, Vη, Vϕ)

Convert a vector `(Vξ, Vη, Vϕ)` in prolate spheroidal components to Cartesian `(Vx, Vy, Vz)`.
"""
function prolate_vector_to_cartesian(a, ξ, η, ϕ, Vξ, Vη, Vϕ)
    hξ, hη, hϕ = scale_factors_prolate(a, ξ, η)
    ex_x, ex_y, ex_z,
    eη_x, eη_y, eη_z,
    eϕ_x, eϕ_y, eϕ_z = prolate_partials(a, ξ, η, ϕ)

    Vx = Vξ*(ex_x/hξ) + Vη*(eη_x/hη) + Vϕ*(eϕ_x/hϕ)
    Vy = Vξ*(ex_y/hξ) + Vη*(eη_y/hη) + Vϕ*(eϕ_y/hϕ)
    Vz = Vξ*(ex_z/hξ) + Vη*(eη_z/hη) + Vϕ*(eϕ_z/hϕ)

    return Vx, Vy, Vz
end

"""
    oblate_vector_to_cartesian(a, ξ, η, ϕ, Vξ, Vη, Vϕ)

Convert a vector `(Vξ, Vη, Vϕ)` in oblate spheroidal components to Cartesian `(Vx, Vy, Vz)`.
"""
function oblate_vector_to_cartesian(a, ξ, η, ϕ, Vξ, Vη, Vϕ)
    hξ, hη, hϕ = scale_factors_oblate(a, ξ, η)

    ex_x, ex_y, ex_z,
    eη_x, eη_y, eη_z,
    eϕ_x, eϕ_y, eϕ_z = oblate_partials(a, ξ, η, ϕ)

    Vx = Vξ*(ex_x/hξ) + Vη*(eη_x/hη) + Vϕ*(eϕ_x/hϕ)
    Vy = Vξ*(ex_y/hξ) + Vη*(eη_y/hη) + Vϕ*(eϕ_y/hϕ)
    Vz = Vξ*(ex_z/hξ) + Vη*(eη_z/hη) + Vϕ*(eϕ_z/hϕ)

    return Vx, Vy, Vz
end


"""
    SpheroidalB{I, T, T2}

Pre-computed data for a single spheroidal wave function with indices `(m, n)` and
parameter `c`. `T2 <: Real` for prolate, `T2 <: Complex` for oblate.

Fields: `m`, `n`, `c`, eigenvalue `λ`, expansion coefficients `dr` and `c2k`.
Constructed via `SpheroidalB(m, n, c)`.
"""
struct SpheroidalB{I, T, T2}
    m::I
    n::I
    c::T2
    λ::T
    dr::Vector{T}
    c2k::Vector{T}
end

function SpheroidalB(m, n, c)
    m > n && throw(ArgumentError("n must be >= m"))
    λ = SpheroidalWaveFunctions.cv_matrix(m, n, c)
    dr = SpheroidalWaveFunctions.compute_dr2_mix(m ,n, c, λ)
    c2k = SpheroidalWaveFunctions.compute_c2k(m, n, dr)
    return SpheroidalB(m, n, c, λ, dr, c2k)
end


# Functions, partial derivative and second partial derivative
function evaluate_angular(b::SpheroidalB{I, T, T}, η) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, ∂S = prolate_angular_ps(m, n, c, λ, c2k, η)
    ∂²S = ((+2η*∂S - (λ - c^2*η^2 + m^2/(1 - η^2))*S) / (1 - η^2))
    return S, ∂S, ∂²S
end

function evaluate_angular(b::SpheroidalB{I, T, Complex{T}}, η) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, ∂S = oblate_angular_ps(m, n, c, λ, c2k, η)
    ∂²S = ((+2η*∂S - (λ + c^2*η^2 - m^2/(1 - η^2))*S) / (1 - η^2))
    return S, ∂S, ∂²S
end

function evaluate_radial1(b::SpheroidalB{I, T, T}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = prolate_radial1(m, n, c, λ, dr, ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 + m^2/(ξ^2 - 1))*R) / (ξ^2 - 1))
    return R, ∂R, ∂²R
end

function evaluate_radial1(b::SpheroidalB{I, T, Complex{T}}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = oblate_radial1(m, n, c, λ, dr, ξ)
    ξ = abs(ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 - m^2/(ξ^2 + 1))*R) / (ξ^2 + 1))
    return R, ∂R, ∂²R
end

function evaluate_radial2(b::SpheroidalB{I, T, T}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = prolate_radial2(m, n, c, λ, dr, ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 + m^2/(ξ^2 - 1))*R) / (ξ^2 - 1))
    return R, ∂R, ∂²R
end

function evaluate_radial2(b::SpheroidalB{I, T, Complex{T}}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = oblate_radial2(m, n, c, λ, dr, ξ)
    ξ = abs(ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 - m^2/(ξ^2 + 1))*R) / (ξ^2 + 1))
    return R, ∂R, ∂²R
end


function evaluate_radial3(b::SpheroidalB{I, T, T}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = prolate_radial1(m, n, c, λ, dr, ξ)
    R2, ∂R2 = prolate_radial2(m, n, c, λ, dr, ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 + m^2/(ξ^2 - 1))*R) / (ξ^2 - 1))
    ∂²R2 = ((-2ξ*∂R2 + (λ - c^2*ξ^2 + m^2/(ξ^2 - 1))*R2) / (ξ^2 - 1))
    return R + im * R2, ∂R + im * ∂R2, ∂²R + im * ∂²R2
end

function evaluate_radial3(b::SpheroidalB{I, T, Complex{T}}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = oblate_radial1(m, n, c, λ, dr, ξ)
    R2, ∂R2 = oblate_radial2(m, n, c, λ, dr, ξ)
    ξ = abs(ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 - m^2/(ξ^2 + 1))*R) / (ξ^2 + 1))
    ∂²R2 = ((-2ξ*∂R2 + (λ - c^2*ξ^2 - m^2/(ξ^2 + 1))*R2) / (ξ^2 + 1))
    return R + im * R2, ∂R + im * ∂R2, ∂²R + im * ∂²R2
end


function evaluate_radial4(b::SpheroidalB{I, T, T}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = prolate_radial1(m, n, c, λ, dr, ξ)
    R2, ∂R2 = prolate_radial2(m, n, c, λ, dr, ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 + m^2/(ξ^2 - 1))*R) / (ξ^2 - 1))
    ∂²R2 = ((-2ξ*∂R2 + (λ - c^2*ξ^2 + m^2/(ξ^2 - 1))*R2) / (ξ^2 - 1))
    return R - im * R2, ∂R - im * ∂R2, ∂²R - im * ∂²R2
end

function evaluate_radial4(b::SpheroidalB{I, T, Complex{T}}, ξ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    R, ∂R = oblate_radial1(m, n, c, λ, dr, ξ)
    R2, ∂R2 = oblate_radial2(m, n, c, λ, dr, ξ)
    ξ = abs(ξ)
    ∂²R = ((-2ξ*∂R + (λ - c^2*ξ^2 - m^2/(ξ^2 + 1))*R) / (ξ^2 + 1))
    ∂²R2 = ((-2ξ*∂R2 + (λ - c^2*ξ^2 - m^2/(ξ^2 + 1))*R2) / (ξ^2 + 1))
    return R - im * R2, ∂R - im * ∂R2, ∂²R - im * ∂²R2
end


abstract type SpheroidalBasis{I, T} end

"""
    ProlateSpheroidalBasis(m_max, n_max, c)

Collection of prolate spheroidal wave functions for all `(m, n)` with
`0 ≤ m ≤ m_max` and `m ≤ n ≤ n_max`. `c` must be real.
"""
struct ProlateSpheroidalBasis{I, T} <: SpheroidalBasis{I, T}
    m_max::I
    n_max::I
    c::T
    basis::Vector{SpheroidalB{I, T, T}}
    function ProlateSpheroidalBasis(m_max, n_max, c::T) where T
        basis = [SpheroidalB(mi, ni, c) for mi in 0:m_max for ni in mi:n_max]
        return new{typeof(m_max), T}(m_max, n_max, c, basis)
    end
end

"""
    OblateSpheroidalBasis(m_max, n_max, c)

Collection of oblate spheroidal wave functions for all `(m, n)` with
`0 ≤ m ≤ m_max` and `m ≤ n ≤ n_max`. `c` must be complex (pass `Complex(c)` for real values).
"""
struct OblateSpheroidalBasis{I, T} <: SpheroidalBasis{I, T}
    m_max::I
    n_max::I
    c::Complex{T}
    basis::Vector{SpheroidalB{I, T, Complex{T}}}
    function OblateSpheroidalBasis(m_max, n_max, c::Complex{T}) where T
        basis = [SpheroidalB(mi, ni, c) for mi in 0:m_max for ni in mi:n_max]
        return new{typeof(m_max), T}(m_max, n_max, c, basis)
    end
end

# Other things

function evaluate(b::SpheroidalB{I, T, T}, η, ϕ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, dS = prolate_angular_ps(m, n, c, λ, c2k, η)
    A = cis(m * ϕ)
    dA = im*m*A
    return S * A, S * dA, dS * A
end

function evaluate(b::SpheroidalB{I, T, Complex{T}}, η, ϕ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, dS = oblate_angular_ps(m, n, c, λ, c2k, η)
    A = cis(m * ϕ)
    dA = im*m*A
    return S * A, S * dA, dS * A
end

function evaluate_real(b::SpheroidalB{I, T, T}, η, ϕ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, dS = prolate_angular_ps(m, n, c, λ, c2k, η)
    sm, cm = sincos(m * ϕ)
    dsm, dcm = m*cm, -m*sm
    return S * sm, S * dsm, dS * sm, S * cm, S * dcm, dS * cm
end

function evaluate_real(b::SpheroidalB{I, T, Complex{T}}, η, ϕ) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, dS =  oblate_angular_ps(m, n, c, λ, c2k, η)
    sm, cm = sincos(m * ϕ)
    dsm, dcm = m*cm, -m*sm
    return S * sm, S * dsm, dS * sm, S * cm, S * dcm, dS * cm
end

"""
    compute_angular_derivatives(basis::SpheroidalBasis, r) -> (S, ∂S, ∂²S)

Batch evaluation of the real angular spheroidal function and its first and second
derivatives with respect to η, for all modes in `basis` and all points in `r`.

`r` is a vector of η values. Returns three real matrices of size `(length(r), N)`,
where `N = length(basis.basis)`.

Use this when the ϕ-dependence is handled separately — in particular, to supply
`S, ∂S, ∂²S` to the N vector functions, which require second derivatives.

See also `compute_angular_with_derivatives` for the full complex mode function
ψ = S(η)·e^{imϕ} with its (η, ϕ) gradient.
"""
function compute_angular_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    N = length(basis.basis)
    S   = Matrix{T}(undef, length(r), N)
    ∂S  = Matrix{T}(undef, length(r), N)
    ∂²S = Matrix{T}(undef, length(r), N)
    for b_idx in eachindex(basis.basis)
        for r_idx in eachindex(r)
            η = r[r_idx]
            Sval, ∂Sval, ∂²Sval = evaluate_angular(basis.basis[b_idx], η)
            S[r_idx, b_idx]   = Sval
            ∂S[r_idx, b_idx]  = ∂Sval
            ∂²S[r_idx, b_idx] = ∂²Sval
        end
    end
    return S, ∂S, ∂²S
end

"""
    compute_angular_with_derivatives(basis::SpheroidalBasis, r) -> (ψ, ∇ψ)

Batch evaluation of the full complex mode function ψ = S(η)·e^{imϕ} and its angular
gradient, for all modes in `basis` and all points in `r`.

`r` is a vector of `(η, ϕ)` tuples. Returns:
- `ψ`: complex matrix of size `(length(r), N)` with `ψ[i,j] = S(η)·e^{imϕ}`.
- `∇ψ`: matrix of `SVector{2, Complex}` with `[∂ψ/∂η, ∂ψ/∂ϕ]` at each point.

Use this for expansion matrices and boundary condition integrals where the full
complex mode function is needed.

See also `compute_angular_derivatives` for the real `S, ∂S, ∂²S` needed
by the N vector functions.
"""
function compute_angular_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    N = length(basis.basis)
    ψ  = Matrix{Complex{T}}(undef, length(r), N)
    ∇ψ = Matrix{SVector{2, Complex{T}}}(undef, length(r), N)
    for b_idx in eachindex(basis.basis)
        for r_idx in eachindex(r)
            η, ϕ = r[r_idx]
            v, ∂v∂η, ∂v∂ϕ = evaluate(basis.basis[b_idx], η, ϕ)
            ψ[r_idx, b_idx] = v
            ∇ψ[r_idx, b_idx] = SVector(∂v∂η, ∂v∂ϕ)
        end
    end
    return ψ, ∇ψ
end

function compute_radial1_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    N = length(basis.basis)
    R   = Matrix{Complex{T}}(undef, length(r), N)
    ∂R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis.basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            Rval, ∂Rval, ∂²Rval = evaluate_radial1(basis.basis[b_idx], ξ)
            R[r_idx, b_idx]   = Rval
            ∂R[r_idx, b_idx]  = ∂Rval
            ∂²R[r_idx, b_idx] = ∂²Rval
        end
    end
    return R, ∂R, ∂²R
end

function compute_radial2_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    N = length(basis.basis)
    R   = Matrix{Complex{T}}(undef, length(r), N)
    ∂R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis.basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            Rval, ∂Rval, ∂²Rval = evaluate_radial2(basis.basis[b_idx], ξ)
            R[r_idx, b_idx]   = Rval
            ∂R[r_idx, b_idx]  = ∂Rval
            ∂²R[r_idx, b_idx] = ∂²Rval
        end
    end
    return R, ∂R, ∂²R
end

function compute_radial3_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    N = length(basis.basis)
    R   = Matrix{Complex{T}}(undef, length(r), N)
    ∂R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis.basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            Rval, ∂Rval, ∂²Rval = evaluate_radial3(basis.basis[b_idx], ξ)
            R[r_idx, b_idx]   = Rval
            ∂R[r_idx, b_idx]  = ∂Rval
            ∂²R[r_idx, b_idx] = ∂²Rval
        end
    end
    return R, ∂R, ∂²R
end

function compute_radial4_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    N = length(basis.basis)
    R   = Matrix{Complex{T}}(undef, length(r), N)
    ∂R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis.basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            Rval, ∂Rval, ∂²Rval = evaluate_radial4(basis.basis[b_idx], ξ)
            R[r_idx, b_idx]   = Rval
            ∂R[r_idx, b_idx]  = ∂Rval
            ∂²R[r_idx, b_idx] = ∂²Rval
        end
    end
    return R, ∂R, ∂²R
end