## Spheroidal coordinates functions 

function obl2cart(a, ξ, η, ϕ)
    ρ = a*sqrt(ξ^2 + 1)*sqrt(1 - η^2)
    x = ρ * cos(ϕ)
    y = ρ * sin(ϕ)
    z = a*ξ*η
    return x,y,z
end

function pro2cart(a, ξ, η, ϕ)
    ρ = a * sqrt(ξ^2 - 1) * sqrt(1 - η^2)
    x = ρ * cos(ϕ)
    y = ρ * sin(ϕ)  
    z = a * ξ * η
    return x, y, z
end

function cart2pro(a, x, y, z)
    ρ = sqrt(x^2 + y^2)
    ϕ = atan(y, x)
    
    r = sqrt(x^2 + y^2 + z^2)
    ξ = 0.5 * (sqrt((r + a)^2 + z^2) + sqrt((r - a)^2 + z^2)) / a
    η = z / (a * ξ)
    
    ξ = max(ξ, 1.0)  
    η = clamp(η, -1.0, 1.0)  # η ∈ [-1, 1]
    
    return ξ, η, ϕ
end

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

function spheroidal_parameter(major_axis, minor_axis)
    a = major_axis / 2
    b = minor_axis / 2
    focal_distance = 2 * sqrt(a^2 - b^2)  
    return focal_distance / 2  
end

function scale_factors_prolate(a, ξ, η)
    hξ = a * sqrt((ξ^2 - 1) * (1 - η^2)) / sqrt(ξ^2 - η^2)
    hη = a * sqrt((ξ^2 - 1) * (1 - η^2)) / sqrt(ξ^2 - η^2)
    hϕ = a * sqrt((ξ^2 - 1) * (1 - η^2))
    return hξ, hη, hϕ
end

function scale_factors_oblate(a, ξ, η)
    hξ = a * sqrt((ξ^2 + 1) * (1 - η^2)) / sqrt(ξ^2 + η^2)
    hη = a * sqrt((ξ^2 + 1) * (1 - η^2)) / sqrt(ξ^2 + η^2)
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


struct SpheroidalB{I, T, T2}
    m::I
    n::I
    c::T2
    λ::T
    dr::Vector{T}
    c2k::Vector{T}
    function SpheroidalB(m, n, c)
        m > n && throw(ArgumentError("n must be >= m"))
        λ = SpheroidalWaveFunctions.cv_matrix(m, n, c)
        dr = SpheroidalWaveFunctions.compute_dr2_mix(m ,n, c, λ)
        c2k = SpheroidalWaveFunctions.compute_c2k(m, n, dr)
        return new{I, T, T2}(m, n, c, λ, dr, c2k)
    end
end


# Functions, partial derivative and second partial derivative
function evaluate_angular(b::SpheroidalB{I, T, T}, η) where {I, T}
    (; m, n, c, λ, dr, c2k) = b
    S, ∂S = prolate_angular_ps(m, n, c, λ, c2k, ξ)
    ∂²S = ((+2η*∂S - (λ - c^2*η^2 + m^2/(ξ^2 - 1))*S) / (1 - η^2))
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
    R, ∂R = oblate_angular_ps(m, n, c, λ, dr, ξ)
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
    R, ∂R = oblate_radial2(m, n, c, λ, ξ)
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

struct ProlateSpheroidalBasis{I, T} <: SpheroidalBasis{I, T}
    m_max::I
    n_max::I
    c::T
    basis::Vector{SpheroidalB{I, T}}
    function SpheroidalBasis(m_max, n_max, c)
        basis = [SpheroidalB(mi, ni, c) for mi in 0:m_max for ni in mi:n_max]
        return new(m_max, n_max, c, basis)
    end
end

struct OblateSpheroidalBasis{I, T} <: SpheroidalBasis{I, T}
    m_max::I
    n_max::I
    c::Complex{T}
    basis::Vector{SpheroidalB{I, T}}
    function SpheroidalBasis(m_max, n_max, c)
        basis = [SpheroidalB(mi, ni, c) for mi in 0:m_max for ni in mi:n_max]
        return new(m_max, n_max, c, basis)
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


function compute_angular_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    m_max = basis.m_max
    n_max = basis.n_max
    m = m_max+1
    n = n_max+1
    N = m*n
    b = basis
    # idx = n_max * (m_idx - 1) + n_idx
    ψ  = Matrix{Complex{T}}(undef, length(r), N)
    ∇ψ = Matrix{SVector{2, Complex{T}}}(undef, length(r), N)
    for b_idx in eachindex(basis)
        for r_idx in eachindex(r)
            η, ϕ = r[r_idx]
            v, ∂v∂η, ∂v∂ϕ  = evaluate(b[idx], η, ϕ)
            ψ[r_idx, b_idx] = v
            ∇ψ[r_idx, b_idx] = SVector(∂v∂η, ∂v∂ϕ)
        end
    end
    return ψ, ∇ψ
end

function compute_radial1_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    m_max = basis.m_max
    n_max = basis.n_max
    m = m_max+1
    n = n_max+1
    N = m*n
    b = basis
    R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂R = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            r, ∂r, ∂²r  = evaluate_radial1(b[idx], ξ)
            R[r_idx, b_idx] = v
            ∂R[r_idx, b_idx] = ∂r
            ∂²R[r_idx, b_idx] = ∂²r
        end
    end
    return R, ∂R, ∂²R
end

function compute_radial2_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    m_max = basis.m_max
    n_max = basis.n_max
    m = m_max+1
    n = n_max+1
    N = m*n
    b = basis
    R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂R = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            r, ∂r, ∂²r  = evaluate_radial2(b[idx], ξ)
            R[r_idx, b_idx] = v
            ∂R[r_idx, b_idx] = ∂r
            ∂²R[r_idx, b_idx] = ∂²r
        end
    end
    return R, ∂R, ∂²R
end

function compute_radial3_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    m_max = basis.m_max
    n_max = basis.n_max
    m = m_max+1
    n = n_max+1
    N = m*n
    b = basis
    R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂R = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            r, ∂r, ∂²r  = evaluate_radial3(b[idx], ξ)
            R[r_idx, b_idx] = v
            ∂R[r_idx, b_idx] = ∂r
            ∂²R[r_idx, b_idx] = ∂²r
        end
    end
    return R, ∂R, ∂²R
end
function compute_radial4_with_derivatives(basis::SpheroidalBasis{I, T}, r) where {I, T}
    m_max = basis.m_max
    n_max = basis.n_max
    m = m_max+1
    n = n_max+1
    N = m*n
    b = basis
    R  = Matrix{Complex{T}}(undef, length(r), N)
    ∂R = Matrix{Complex{T}}(undef, length(r), N)
    ∂²R = Matrix{Complex{T}}(undef, length(r), N)
    for b_idx in eachindex(basis)
        for r_idx in eachindex(r)
            ξ = r[r_idx]
            r, ∂r, ∂²r  = evaluate_radial4(b[idx], ξ)
            R[r_idx, b_idx] = v
            ∂R[r_idx, b_idx] = ∂r
            ∂²R[r_idx, b_idx] = ∂²r
        end
    end
    return R, ∂R, ∂²R
end