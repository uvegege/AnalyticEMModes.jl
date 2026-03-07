include("./spheroidal_utils.jl")
include("./spheroidal_vectors.jl")


kc_spheroidal() = 0.0

spheroidal_families() = (:x, :y, :z, :r)

@inline _is_oblate_basis(::SpheroidalB{I, T, T}) where {I, T} = false
@inline _is_oblate_basis(::SpheroidalB{I, T, Complex{T}}) where {I, T} = true

function _eval_M(family::Symbol, ξ, η, ϕ, S, dS, R, dR, m, d, even, sph_type)
    if family == :x
        return Mˣ(ξ, η, ϕ, S, dS, R, dR, m, d, even, sph_type)
    elseif family == :y
        return Mʸ(ξ, η, ϕ, S, dS, R, dR, m, d, even, sph_type)
    elseif family == :z
        return Mᶻ(ξ, η, ϕ, S, dS, R, dR, m, d, even, sph_type)
    else
        return Mʳ(ξ, η, ϕ, S, dS, R, dR, m, d, even, sph_type)
    end
end

function _eval_N(family::Symbol, ξ, η, ϕ, S, dS, d²S, R, dR, d²R, m, d, k, even, oblate)
    if family == :x
        return Nˣ(ξ, η, ϕ, S, dS, d²S, R, dR, d²R, m, d, k, even, oblate)
    elseif family == :y
        return Nʸ(ξ, η, ϕ, S, dS, d²S, R, dR, d²R, m, d, k, even, oblate)
    elseif family == :z
        return Nᶻ(ξ, η, ϕ, S, dS, d²S, R, dR, d²R, m, d, k, even, oblate)
    else
        return Nʳ(ξ, η, ϕ, S, dS, d²S, R, dR, d²R, m, d, k, even, oblate)
    end
end

function _evaluate_radial(basis, ξ, radial::Int)
    if radial == 1
        return evaluate_radial1(basis, ξ)
    elseif radial == 2
        return evaluate_radial2(basis, ξ)
    elseif radial == 3
        return evaluate_radial3(basis, ξ)
    else
        return evaluate_radial4(basis, ξ)
    end
end

function _compute_radial(basis::SpheroidalBasis{I, T}, ξs, radial::Int) where {I, T}
    if radial == 1
        return compute_radial1_with_derivatives(basis, ξs)
    elseif radial == 2
        return compute_radial2_with_derivatives(basis, ξs)
    elseif radial == 3
        return compute_radial3_with_derivatives(basis, ξs)
    else
        return compute_radial4_with_derivatives(basis, ξs)
    end
end

"""
    mn_spheroidal_vector(ξ, η, ϕ, mode::SpheroidalB, k;
                         family=:z, even=true, oblate=_is_oblate_basis(mode), radial=4)

Evaluate one spheroidal mode and return `(Mξ, Mη, Mϕ, Nξ, Nη, Nϕ)` in local spheroidal components.

# Keyword arguments
- `family`: vector family — `:x`, `:y`, `:z`, or `:r`. Default: `:z`.
- `even`: parity — `true` for cos(mϕ) (even), `false` for sin(mϕ) (odd). Default: `true`.
- `oblate`: use oblate conventions if `true`, prolate if `false`. Inferred from `mode` type by default.
- `radial`: radial function type (1=R₁, 2=R₂, 3=R₃ incoming, 4=R₄ outgoing). Default: `4`.
"""
function mn_spheroidal_vector(ξ, η, ϕ, mode::SpheroidalB, k;
                              family::Symbol = :z, even::Bool = true,
                              oblate::Bool = _is_oblate_basis(mode),
                              radial::Int = 4)
    S, dS, d²S = evaluate_angular(mode, η)
    R, dR, d²R = _evaluate_radial(mode, ξ, radial)
    d = mode.c / k
    sph_type = oblate ? :oblate : :prolate
    Mξ, Mη, Mϕ = _eval_M(family, ξ, η, ϕ, S, dS, R, dR, mode.m, d, even, sph_type)
    Nξ, Nη, Nϕ = _eval_N(family, ξ, η, ϕ, S, dS, d²S, R, dR, d²R, mode.m, d, k, even, oblate)
    return (Mξ, Mη, Mϕ, Nξ, Nη, Nϕ)
end

function mn_spheroidal_vector(p::SVector{3, T}, mode::SpheroidalB, k; kwargs...) where {T}
    return mn_spheroidal_vector(p[1], p[2], p[3], mode, k; kwargs...)
end

"""
    m_spheroidal_vector(ξ, η, ϕ, mode::SpheroidalB, k; kwargs...)

Evaluate one spheroidal mode and return `(Mξ, Mη, Mϕ)`.
"""
function m_spheroidal_vector(ξ, η, ϕ, mode::SpheroidalB, k; kwargs...)
    Mξ, Mη, Mϕ, _, _, _ = mn_spheroidal_vector(ξ, η, ϕ, mode, k; kwargs...)
    return (Mξ, Mη, Mϕ)
end

function m_spheroidal_vector(p::SVector{3, T}, mode::SpheroidalB, k; kwargs...) where {T}
    return m_spheroidal_vector(p[1], p[2], p[3], mode, k; kwargs...)
end

"""
    n_spheroidal_vector(ξ, η, ϕ, mode::SpheroidalB, k; kwargs...)

Evaluate one spheroidal mode and return `(Nξ, Nη, Nϕ)`.
"""
function n_spheroidal_vector(ξ, η, ϕ, mode::SpheroidalB, k; kwargs...)
    _, _, _, Nξ, Nη, Nϕ = mn_spheroidal_vector(ξ, η, ϕ, mode, k; kwargs...)
    return (Nξ, Nη, Nϕ)
end

function n_spheroidal_vector(p::SVector{3, T}, mode::SpheroidalB, k; kwargs...) where {T}
    return n_spheroidal_vector(p[1], p[2], p[3], mode, k; kwargs...)
end

"""
    mn_spheroidal_vectors(points, basis::SpheroidalBasis, k;
                          family=:z, even=true, oblate=nothing, radial=4)

Batch evaluation of M and N spheroidal vectors over all points and all modes in `basis`.
Returns a matrix `A[p, i] = (Mξ, Mη, Mϕ, Nξ, Nη, Nϕ)` in local spheroidal components.

`points` can be a vector of `SVector{3}` or any indexable collection with `p[1]=ξ, p[2]=η, p[3]=ϕ`.

# Keyword arguments
- `family`: vector family — `:x`, `:y`, `:z`, or `:r`. Default: `:z`.
- `even`: parity — `true` for cos(mϕ), `false` for sin(mϕ). Default: `true`.
- `oblate`: override oblate/prolate selection (`nothing` infers from basis type).
- `radial`: radial function type (1=R₁, 2=R₂, 3=R₃ incoming, 4=R₄ outgoing). Default: `4`.
"""
function mn_spheroidal_vectors(points, basis::SpheroidalBasis{I, T}, k;
                               family::Symbol = :z, even::Bool = true,
                               oblate = nothing, radial::Int = 4) where {I, T}
    npts = length(points)
    nmodes = length(basis.basis)
    A = Matrix{NTuple{6, Complex{T}}}(undef, npts, nmodes)
    _mn_spheroidal_vectors!(A, points, basis, k; family = family, even = even, oblate = oblate, radial = radial)
    return A
end

function mn_spheroidal_vectors(points::AbstractVector{<:SVector{3}}, basis::SpheroidalBasis{I, T}, k;
                               family::Symbol = :z, even::Bool = true,
                               oblate = nothing, radial::Int = 4) where {I, T}
    npts = length(points)
    nmodes = length(basis.basis)
    A = Matrix{NTuple{6, Complex{T}}}(undef, npts, nmodes)
    _mn_spheroidal_vectors!(A, points, basis, k; family = family, even = even, oblate = oblate, radial = radial)
    return A
end

function _mn_spheroidal_vectors!(A, points, basis::SpheroidalBasis{I, T}, k;
                                 family::Symbol = :z, even::Bool = true,
                                 oblate = nothing, radial::Int = 4) where {I, T}
    sph_oblate = isnothing(oblate) ? (basis isa OblateSpheroidalBasis) : oblate
    sph_type = sph_oblate ? :oblate : :prolate

    ηϕ = [(p[2], p[3]) for p in points]
    ξs = [p[1] for p in points]
    ηs = [p[2] for p in points]
    ψ, ∇ψ = compute_angular_with_derivatives(basis, ηϕ)
    Sm, dSm, d²Sm = compute_angular_derivatives(basis, ηs)
    Rm, dRm, d²Rm = _compute_radial(basis, ξs, radial)

    for midx in eachindex(basis.basis)
        b = basis.basis[midx]
        d = b.c / k
        for pidx in eachindex(points)
            ξ, η, ϕ = points[pidx]
            S = ψ[pidx, midx]
            dS = ∇ψ[pidx, midx][1]
            R   = Rm[pidx, midx]
            dR  = dRm[pidx, midx]
            d²R = d²Rm[pidx, midx]
            S_r, dS_r, d²S_r = Sm[pidx, midx], dSm[pidx, midx], d²Sm[pidx, midx]
            Mξ, Mη, Mϕ = _eval_M(family, ξ, η, ϕ, S, dS, R, dR, b.m, d, even, sph_type)
            Nξ, Nη, Nϕ = _eval_N(family, ξ, η, ϕ, S_r, dS_r, d²S_r, R, dR, d²R, b.m, d, k, even, sph_oblate)
            A[pidx, midx] = (Mξ, Mη, Mϕ, Nξ, Nη, Nϕ)
        end
    end
    return nothing
end

"""
    m_spheroidal_vectors(points, basis::SpheroidalBasis, k; kwargs...)

Batch evaluation returning only `(Mξ, Mη, Mϕ)`.
"""
function m_spheroidal_vectors(points, basis::SpheroidalBasis{I, T}, k; kwargs...) where {I, T}
    npts = length(points)
    nmodes = length(basis.basis)
    A = Matrix{NTuple{3, Complex{T}}}(undef, npts, nmodes)
    B = mn_spheroidal_vectors(points, basis, k; kwargs...)
    for idx in eachindex(B)
        A[idx] = (B[idx][1], B[idx][2], B[idx][3])
    end
    return A
end

"""
    n_spheroidal_vectors(points, basis::SpheroidalBasis, k; kwargs...)

Batch evaluation returning only `(Nξ, Nη, Nϕ)`.
"""
function n_spheroidal_vectors(points, basis::SpheroidalBasis{I, T}, k; kwargs...) where {I, T}
    npts = length(points)
    nmodes = length(basis.basis)
    A = Matrix{NTuple{3, Complex{T}}}(undef, npts, nmodes)
    B = mn_spheroidal_vectors(points, basis, k; kwargs...)
    for idx in eachindex(B)
        A[idx] = (B[idx][4], B[idx][5], B[idx][6])
    end
    return A
end

"""
    mn_spheroidal_vectors(xi, eta, phi, basis::SpheroidalBasis, k; kwargs...)

Convenience overload accepting separate coordinate arrays `xi`, `eta`, `phi`.
See [`mn_spheroidal_vectors`](@ref) for keyword arguments.
"""
function mn_spheroidal_vectors(xi::AbstractArray, eta::AbstractArray, phi::AbstractArray,
                               basis::SpheroidalBasis, k; kwargs...)
    points = to_svector(xi, eta, phi)
    return mn_spheroidal_vectors(points, basis, k; kwargs...)
end

"""
    mn_spheroidal_vectors_mnmax(points, m_max, n_max, c, k;
                                oblate=false, family=:z, even=true, radial=4)

Build a spheroidal basis with all `(m, n)` for `0 ≤ m ≤ m_max`, `m ≤ n ≤ n_max` and
evaluate all M/N vectors in one call. `c` is the spheroidal parameter.

See [`mn_spheroidal_vectors`](@ref) for keyword arguments.
"""
function mn_spheroidal_vectors_mnmax(points, m_max, n_max, c, k;
                                     oblate::Bool = false, family::Symbol = :z,
                                     even::Bool = true, radial::Int = 4)
    basis = oblate ? OblateSpheroidalBasis(m_max, n_max, Complex(c)) :
                     ProlateSpheroidalBasis(m_max, n_max, c)
    return mn_spheroidal_vectors(points, basis, k; family = family, even = even, oblate = oblate, radial = radial)
end
