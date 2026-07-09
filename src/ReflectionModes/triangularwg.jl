"""
    ReflectionTriangle

Abstract type for triangular waveguides whose scalar Helmholtz modes are built
from reflection symmetries.
"""
abstract type ReflectionTriangle end

"""
    EquilateralTriangle(side)

Equilateral triangular cross-section with side length `side`.

The reference triangle has vertices `(0, 0)`, `(side, 0)`, and
`(side / 2, sqrt(3) * side / 2)`.
"""
struct EquilateralTriangle{T} <: ReflectionTriangle
    side::T
end

"""
    RightIsoscelesTriangle(side)

Right isosceles triangular cross-section obtained by cutting a square of side
`side` along the diagonal `x = y`.
"""
struct RightIsoscelesTriangle{T} <: ReflectionTriangle
    side::T
end

"""
    HalfEquilateralTriangle(side)

30-60-90 triangular cross-section represented as half of an equilateral triangle
of side length `side`.
"""
struct HalfEquilateralTriangle{T} <: ReflectionTriangle
    side::T
end

"""
    ReflectionMode(triangle, m, n, kind)

Scalar reflection mode on a triangular waveguide. `kind` must be `:TE` or `:TM`.
"""
struct ReflectionMode{G <: ReflectionTriangle}
    triangle::G
    m::Int
    n::Int
    kind::Symbol
    symmetry::Symbol
end

TE(triangle::ReflectionTriangle, m, n; symmetry=:S) = ReflectionMode(triangle, Int(m), Int(n), :TE, symmetry)

TM(triangle::ReflectionTriangle, m, n; symmetry=:S) = ReflectionMode(triangle, Int(m), Int(n), :TM, symmetry)

ReflectionMode(triangle::ReflectionTriangle, m, n, kind::Symbol) = ReflectionMode(triangle, Int(m), Int(n), kind, :S)

function check_modekind(kind)
    if !(kind == :TE || kind == :TM)
        throw(ArgumentError("mode kind must be :TE or :TM"))
    end
    return nothing
end

function check_symmetry(symmetry)
    if !(symmetry == :S || symmetry == :A)
        throw(ArgumentError("mode symmetry must be :S or :A"))
    end
    return nothing
end

"""
    kc_equilateral(side, m, n)

Cutoff wavenumber of the `(m, n)` equilateral-triangle mode.
"""
kc_equilateral(side, m, n) = 4 * π / (3 * side) * sqrt(m^2 + m * n + n^2)

kc(mode::ReflectionMode{<:EquilateralTriangle}) = kc_equilateral(mode.triangle.side, mode.m, mode.n)

function equilateral_tm_modal_f(u, v, side, m, n)
    return equilateral_tm_modal_f(u, v, side, m, n, :S)
end

function equilateral_tm_modal_f(u, v, side, m, n, symmetry)
    check_symmetry(symmetry)
    h = sqrt(3) * side / 2
    X = u - side / 2
    Y = v - h
    ky1 = (-m - n) * π / h
    ky2 = m * π / h
    ky3 = n * π / h
    kx1 = (m - n) * π / (sqrt(3) * h)
    kx2 = (2 * n + m) * π / (sqrt(3) * h)
    kx3 = (-2 * m - n) * π / (sqrt(3) * h)

    sy1, cy1 = sincos(ky1 * Y)
    sy2, cy2 = sincos(ky2 * Y)
    sy3, cy3 = sincos(ky3 * Y)
    sx1, cx1 = sincos(kx1 * X)
    sx2, cx2 = sincos(kx2 * X)
    sx3, cx3 = sincos(kx3 * X)

    if symmetry == :S
        ψₖ = sy1 * cx1 + sy2 * cx2 + sy3 * cx3
        ∂ψᵢ = -kx1 * sy1 * sx1 - kx2 * sy2 * sx2 - kx3 * sy3 * sx3
        ∂ψⱼ = ky1 * cy1 * cx1 + ky2 * cy2 * cx2 + ky3 * cy3 * cx3
    else
        ψₖ = sy1 * sx1 + sy2 * sx2 + sy3 * sx3
        ∂ψᵢ = kx1 * sy1 * cx1 + kx2 * sy2 * cx2 + kx3 * sy3 * cx3
        ∂ψⱼ = ky1 * cy1 * sx1 + ky2 * cy2 * sx2 + ky3 * cy3 * sx3
    end
    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end

function equilateral_te_modal_f(u, v, side, m, n)
    return equilateral_te_modal_f(u, v, side, m, n, :S)
end

function equilateral_te_modal_f(u, v, side, m, n, symmetry)
    check_symmetry(symmetry)
    h = sqrt(3) * side / 2
    X = u - side / 2
    Y = v - h
    ky1 = (-m - n) * π / h
    ky2 = m * π / h
    ky3 = n * π / h
    kx1 = (m - n) * π / (sqrt(3) * h)
    kx2 = (2 * n + m) * π / (sqrt(3) * h)
    kx3 = (-2 * m - n) * π / (sqrt(3) * h)

    sy1, cy1 = sincos(ky1 * Y)
    sy2, cy2 = sincos(ky2 * Y)
    sy3, cy3 = sincos(ky3 * Y)
    sx1, cx1 = sincos(kx1 * X)
    sx2, cx2 = sincos(kx2 * X)
    sx3, cx3 = sincos(kx3 * X)

    if symmetry == :S
        ψₖ = cy1 * cx1 + cy2 * cx2 + cy3 * cx3
        ∂ψᵢ = -kx1 * cy1 * sx1 - kx2 * cy2 * sx2 - kx3 * cy3 * sx3
        ∂ψⱼ = -ky1 * sy1 * cx1 - ky2 * sy2 * cx2 - ky3 * sy3 * cx3
    else
        ψₖ = cy1 * sx1 + cy2 * sx2 + cy3 * sx3
        ∂ψᵢ = kx1 * cy1 * cx1 + kx2 * cy2 * cx2 + kx3 * cy3 * cx3
        ∂ψⱼ = -ky1 * sy1 * sx1 - ky2 * sy2 * sx2 - ky3 * sy3 * sx3
    end
    return (∂ψᵢ, ∂ψⱼ, ψₖ)
end

function equilateral_modal_f(x, y, triangle::EquilateralTriangle, m, n, kind)
    return equilateral_modal_f(x, y, triangle, m, n, kind, :S)
end

function equilateral_modal_f(x, y, triangle::EquilateralTriangle, m, n, kind, symmetry)
    check_modekind(kind)
    check_symmetry(symmetry)
    return kind == :TE ? equilateral_te_modal_f(x, y, triangle.side, m, n, symmetry) : equilateral_tm_modal_f(x, y, triangle.side, m, n, symmetry)
end

equilateral_modal_f(x, y, side, m, n, kind) = equilateral_modal_f(x, y, EquilateralTriangle(side), m, n, kind)
equilateral_modal_f(x, y, side, m, n, kind, symmetry) = equilateral_modal_f(x, y, EquilateralTriangle(side), m, n, kind, symmetry)

ψ(mode::ReflectionMode{<:EquilateralTriangle}, x, y) = equilateral_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind, mode.symmetry)[3]

∇ψ(mode::ReflectionMode{<:EquilateralTriangle}, x, y) = equilateral_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind, mode.symmetry)[1:2]

"""
    kc_right_isosceles(side, m, n)

Cutoff wavenumber of the `(m, n)` right-isosceles-triangle mode.
"""
kc_right_isosceles(side, m, n) = π / side * hypot(m, n)

kc(mode::ReflectionMode{<:RightIsoscelesTriangle}) = kc_right_isosceles(mode.triangle.side, mode.m, mode.n)

function square_tm_f(u, v, side, m, n)
    ∂x, ∂y, ψₖ = tm_rwg_modal_f(u, v, side, side, m, n)
    return (∂x, ∂y, ψₖ)
end

function square_te_f(u, v, side, m, n)
    ∂x, ∂y, ψₖ = te_rwg_modal_f(u, v, side, side, m, n)
    return (∂x, ∂y, ψₖ)
end

function right_isosceles_modal_f(x, y, triangle::RightIsoscelesTriangle, m, n, kind)
    check_modekind(kind)

    if kind == :TE
        ∂mn_x, ∂mn_y, ψmn = square_te_f(x, y, triangle.side, m, n)
        ∂nm_x, ∂nm_y, ψnm = square_te_f(x, y, triangle.side, n, m)
        s = one(ψmn)
    else
        ∂mn_x, ∂mn_y, ψmn = square_tm_f(x, y, triangle.side, m, n)
        ∂nm_x, ∂nm_y, ψnm = square_tm_f(x, y, triangle.side, n, m)
        s = -one(ψmn)
    end

    ∂x = ∂mn_x + s * ∂nm_x
    ∂y = ∂mn_y + s * ∂nm_y
    ψₖ = ψmn + s * ψnm
    return (∂x, ∂y, ψₖ)
end

right_isosceles_modal_f(x, y, side, m, n, kind) = right_isosceles_modal_f(x, y, RightIsoscelesTriangle(side), m, n, kind)

ψ(mode::ReflectionMode{<:RightIsoscelesTriangle}, x, y) = right_isosceles_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind)[3]

∇ψ(mode::ReflectionMode{<:RightIsoscelesTriangle}, x, y) = right_isosceles_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind)[1:2]

"""
    kc_half_equilateral(side, m, n)

Cutoff wavenumber inherited from the corresponding equilateral mode.
"""
kc_half_equilateral(side, m, n) = kc_equilateral(side, m, n)

kc(mode::ReflectionMode{<:HalfEquilateralTriangle}) = kc_half_equilateral(mode.triangle.side, mode.m, mode.n)

function scalar_integral(mode::ReflectionMode{<:EquilateralTriangle})
    area = sqrt(3) * mode.triangle.side^2 / 4   
    mode.m + mode.n == 0 && return zero(area)
    mode.symmetry == :A && mode.m == mode.n && return zero(area)
    mode.kind == :TM && (mode.m == 0 || mode.n == 0) && return zero(area)
    return mode.m == mode.n || mode.m == 0 || mode.n == 0 ? 3 * area / 2 : 3 * area / 4
end

function scalar_integral(mode::ReflectionMode{<:RightIsoscelesTriangle})
    area = mode.triangle.side^2 / 2
    mode.m + mode.n == 0 && return zero(area)
    mode.kind == :TM && (mode.m == 0 || mode.n == 0 || mode.m == mode.n) && return zero(area)
    return mode.m == mode.n || mode.m == 0 || mode.n == 0 ? area : area / 2
end

function scalar_integral(mode::ReflectionMode{<:HalfEquilateralTriangle})
    area = sqrt(3) * mode.triangle.side^2 / 8
    mode.m + mode.n == 0 && return zero(area)
    mode.kind == :TM && (mode.m == 0 || mode.n == 0 || mode.m == mode.n) && return zero(area)
    return mode.m == mode.n || mode.m == 0 || mode.n == 0 ? 3 * area / 2 : 3 * area / 4
end

function triangular_modal_power(mode::ReflectionMode, k, β, f, μᵣ, εᵣ)
    ω = 2 * π * f
    factor = mode.kind == :TE ? ω * μᵣ * _μₒ * β / k^2 : ω * εᵣ * _εₒ * β / k^2
    return 0.5 * factor * scalar_integral(mode)
end

te_normalization(mode::ReflectionMode, β, f, μᵣ, εᵣ) = sqrt(1 / triangular_modal_power(mode, kc(mode), β, f, μᵣ, εᵣ))

tm_normalization(mode::ReflectionMode, β, f, μᵣ, εᵣ) = sqrt(1 / triangular_modal_power(mode, kc(mode), β, f, μᵣ, εᵣ))

"""
    te_normalization_equilateral(side, m, n, symmetry, kc, β, f, μᵣ, εᵣ)

Normalization factor for equilateral-triangle TE modes to achieve unit power.
"""
function te_normalization_equilateral(side, m, n, symmetry, kc, β, f, μᵣ, εᵣ)
    mode = ReflectionMode(EquilateralTriangle(side), m, n, :TE, symmetry)
    return sqrt(1 / triangular_modal_power(mode, kc, β, f, μᵣ, εᵣ))
end

te_normalization_equilateral(side, m, n, kc, β, f, μᵣ, εᵣ) = te_normalization_equilateral(side, m, n, :S, kc, β, f, μᵣ, εᵣ)

"""
    tm_normalization_equilateral(side, m, n, symmetry, kc, β, f, μᵣ, εᵣ)

Normalization factor for equilateral-triangle TM modes to achieve unit power.
"""
function tm_normalization_equilateral(side, m, n, symmetry, kc, β, f, μᵣ, εᵣ)
    mode = ReflectionMode(EquilateralTriangle(side), m, n, :TM, symmetry)
    return sqrt(1 / triangular_modal_power(mode, kc, β, f, μᵣ, εᵣ))
end

tm_normalization_equilateral(side, m, n, kc, β, f, μᵣ, εᵣ) = tm_normalization_equilateral(side, m, n, :S, kc, β, f, μᵣ, εᵣ)

"""
    te_normalization_right_isosceles(side, m, n, kc, β, f, μᵣ, εᵣ)

Normalization factor for right-isosceles-triangle TE modes to achieve unit power.
"""
function te_normalization_right_isosceles(side, m, n, kc, β, f, μᵣ, εᵣ)
    mode = ReflectionMode(RightIsoscelesTriangle(side), m, n, :TE)
    return sqrt(1 / triangular_modal_power(mode, kc, β, f, μᵣ, εᵣ))
end

"""
    tm_normalization_right_isosceles(side, m, n, kc, β, f, μᵣ, εᵣ)

Normalization factor for right-isosceles-triangle TM modes to achieve unit power.
"""
function tm_normalization_right_isosceles(side, m, n, kc, β, f, μᵣ, εᵣ)
    mode = ReflectionMode(RightIsoscelesTriangle(side), m, n, :TM)
    return sqrt(1 / triangular_modal_power(mode, kc, β, f, μᵣ, εᵣ))
end

"""
    te_normalization_half_equilateral(side, m, n, kc, β, f, μᵣ, εᵣ)

Normalization factor for half-equilateral-triangle TE modes to achieve unit power.
"""
function te_normalization_half_equilateral(side, m, n, kc, β, f, μᵣ, εᵣ)
    mode = ReflectionMode(HalfEquilateralTriangle(side), m, n, :TE)
    return sqrt(1 / triangular_modal_power(mode, kc, β, f, μᵣ, εᵣ))
end

"""
    tm_normalization_half_equilateral(side, m, n, kc, β, f, μᵣ, εᵣ)

Normalization factor for half-equilateral-triangle TM modes to achieve unit power.
"""
function tm_normalization_half_equilateral(side, m, n, kc, β, f, μᵣ, εᵣ)
    mode = ReflectionMode(HalfEquilateralTriangle(side), m, n, :TM)
    return sqrt(1 / triangular_modal_power(mode, kc, β, f, μᵣ, εᵣ))
end

"""
    first_n_modes_equilateral(N, side)

Returns the first `N` distinct scalar modes for an equilateral triangular
waveguide. Permuted index pairs that produce the same scalar pattern are
represented only once.

The returned tuples are `(kind, m, n, symmetry, kc)`.
"""
function first_n_modes_equilateral(N, side)
    N <= 0 && return Tuple{Symbol, Int, Int, Symbol, Float64}[]
    modes = Tuple{Symbol, Int, Int, Symbol, Float64}[]
    max_index = 1
    while true
        empty!(modes)
        # The equilateral symmetry group has no orientation family larger than
        # two dimensions; `:S` and `:A` span the generic pair, while `(n, m)`
        # repeats `(m, n)` up to sign.
        for m in 0:max_index, n in m:max_index
            m + n == 0 && continue
            push!(modes, (:TE, m, n, :S, kc_equilateral(side, m, n)))
            m < n && push!(modes, (:TE, m, n, :A, kc_equilateral(side, m, n)))
            if m >= 1
                push!(modes, (:TM, m, n, :S, kc_equilateral(side, m, n)))
                m < n && push!(modes, (:TM, m, n, :A, kc_equilateral(side, m, n)))
            end
        end
        sort!(modes, by=x -> (x[5], x[1] == :TE ? 0 : 1, x[4] == :S ? 0 : 1, x[2], x[3]))
        length(modes) >= N && modes[N][5] < kc_equilateral(side, 0, max_index + 1) && break
        max_index *= 2
    end
    return modes[1:N]
end

"""
    first_n_modes_right_isosceles(N, side)

Returns the first `N` distinct scalar modes for a right isosceles triangular
waveguide. Symmetric TE and antisymmetric TM square modes are listed once per
distinct triangle mode.

The returned tuples are `(kind, m, n, kc)`.
"""
function first_n_modes_right_isosceles(N, side)
    N <= 0 && return Tuple{Symbol, Int, Int, Float64}[]
    modes = Tuple{Symbol, Int, Int, Float64}[]
    max_index = 1
    while true
        empty!(modes)
        for m in 0:max_index, n in m:max_index
            m + n > 0 && push!(modes, (:TE, m, n, kc_right_isosceles(side, m, n)))
            m >= 1 && m < n && push!(modes, (:TM, m, n, kc_right_isosceles(side, m, n)))
        end
        sort!(modes, by=x -> (x[4], x[1] == :TE ? 0 : 1, x[2], x[3]))
        length(modes) >= N && modes[N][4] < kc_right_isosceles(side, 0, max_index + 1) && break
        max_index *= 2
    end
    return modes[1:N]
end

"""
    first_n_modes_half_equilateral(N, side)

Returns the first `N` distinct scalar modes for a half-equilateral triangular
waveguide. The TE modes come from symmetric equilateral modes and the TM modes
from antisymmetric equilateral modes.

The returned tuples are `(kind, m, n, kc)`.
"""
function first_n_modes_half_equilateral(N, side)
    N <= 0 && return Tuple{Symbol, Int, Int, Float64}[]
    modes = Tuple{Symbol, Int, Int, Float64}[]
    max_index = 1
    while true
        empty!(modes)
        for m in 0:max_index, n in m:max_index
            m + n > 0 && push!(modes, (:TE, m, n, kc_half_equilateral(side, m, n)))
            m >= 1 && m < n && push!(modes, (:TM, m, n, kc_half_equilateral(side, m, n)))
        end
        sort!(modes, by=x -> (x[4], x[1] == :TE ? 0 : 1, x[2], x[3]))
        length(modes) >= N && modes[N][4] < kc_half_equilateral(side, 0, max_index + 1) && break
        max_index *= 2
    end
    return modes[1:N]
end

function half_equilateral_modal_f(x, y, triangle::HalfEquilateralTriangle, m, n, kind)
    check_modekind(kind)
    if kind == :TE
        ∂x, ∂y, ψₖ = equilateral_te_modal_f(x + triangle.side / 2, y, triangle.side, m, n, :S)
    else
        ∂x, ∂y, ψₖ = equilateral_tm_modal_f(x + triangle.side / 2, y, triangle.side, m, n, :A)
    end

    return (∂x, ∂y, ψₖ)
end

half_equilateral_modal_f(x, y, side, m, n, kind) = half_equilateral_modal_f(x, y, HalfEquilateralTriangle(side), m, n, kind)

ψ(mode::ReflectionMode{<:HalfEquilateralTriangle}, x, y) = half_equilateral_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind)[3]

∇ψ(mode::ReflectionMode{<:HalfEquilateralTriangle}, x, y) = half_equilateral_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind)[1:2]

function te_reflection_fields_from_modal(∂ψᵢ, ∂ψⱼ, ψₖ, c_e, c_h)
    Ex = -c_e * ∂ψᵢ
    Ey = +c_e * ∂ψⱼ
    Ez = zero(Ex)
    Hx = -c_h * ∂ψⱼ
    Hy = -c_h * ∂ψᵢ
    Hz = -im * ψₖ
    return (Ex, Ey, Ez, Hx, Hy, Hz)
end

function tm_reflection_fields_from_modal(∂ψᵢ, ∂ψⱼ, ψₖ, c_e, c_h)
    Ex = -c_e * ∂ψᵢ
    Ey = -c_e * ∂ψⱼ
    Ez = -im * ψₖ
    Hx = +c_h * ∂ψⱼ
    Hy = -c_h * ∂ψᵢ
    Hz = zero(Ex)
    return (Ex, Ey, Ez, Hx, Hy, Hz)
end

function te_reflection_fields(mode::ReflectionMode, x, y, c_e, c_h)
    mode.kind == :TE || throw(ArgumentError("mode kind must be :TE"))
    ∂ψᵢ, ∂ψⱼ, ψₖ = modal_f(mode, x, y)
    return te_reflection_fields_from_modal(∂ψᵢ, ∂ψⱼ, ψₖ, c_e, c_h)
end

function te_reflection_fields(mode::ReflectionMode, x, y, f, μᵣ, εᵣ)
    k = kc(mode)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(k, β, f, μᵣ, εᵣ)
    return te_reflection_fields(mode, x, y, c_e, c_h)
end

function tm_reflection_fields(mode::ReflectionMode, x, y, c_e, c_h)
    mode.kind == :TM || throw(ArgumentError("mode kind must be :TM"))
    ∂ψᵢ, ∂ψⱼ, ψₖ = modal_f(mode, x, y)
    return tm_reflection_fields_from_modal(∂ψᵢ, ∂ψⱼ, ψₖ, c_e, c_h)
end

function tm_reflection_fields(mode::ReflectionMode, x, y, f, μᵣ, εᵣ)
    k = kc(mode)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(k, β, f, μᵣ, εᵣ)
    return tm_reflection_fields(mode, x, y, c_e, c_h)
end

function modal_f(mode::ReflectionMode{<:EquilateralTriangle}, x, y)
    return equilateral_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind, mode.symmetry)
end

function modal_f(mode::ReflectionMode{<:RightIsoscelesTriangle}, x, y)
    return right_isosceles_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind)
end

function modal_f(mode::ReflectionMode{<:HalfEquilateralTriangle}, x, y)
    return half_equilateral_modal_f(x, y, mode.triangle, mode.m, mode.n, mode.kind)
end

function te_reflection_fields(mode::ReflectionMode, x::AbstractArray{T, N}, y::AbstractArray{T, N}, f, μᵣ, εᵣ) where {T, N}
    k = kc(mode)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(k, β, f, μᵣ, εᵣ)
    fields = similar(x, NTuple{6, Complex{T}})
    for idx in eachindex(x)
        fields[idx] = te_reflection_fields(mode, x[idx], y[idx], c_e, c_h)
    end
    return fields
end

function tm_reflection_fields(mode::ReflectionMode, x::AbstractArray{T, N}, y::AbstractArray{T, N}, f, μᵣ, εᵣ) where {T, N}
    k = kc(mode)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(k, β, f, μᵣ, εᵣ)
    fields = similar(x, NTuple{6, Complex{T}})
    for idx in eachindex(x)
        fields[idx] = tm_reflection_fields(mode, x[idx], y[idx], c_e, c_h)
    end
    return fields
end

function array_fields(fieldfun::F, x::AbstractArray{T, N}, y::AbstractArray{T, N}, side, m, n, f, μᵣ, εᵣ) where {F, T, N}
    fields = similar(x, NTuple{6, Complex{T}})
    for idx in eachindex(x)
        fields[idx] = fieldfun(x[idx], y[idx], side, m, n, f, μᵣ, εᵣ)
    end
    return fields
end

"""
    te_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)

Electric and magnetic field components of an equilateral-triangle TE mode.
Returns `(Ex, Ey, Ez, Hx, Hy, Hz)`.
"""
te_equilateral_fields(x, y, side, m, n, c_e, c_h) = te_reflection_fields_from_modal(equilateral_modal_f(x, y, side, m, n, :TE, :S)..., c_e, c_h)

function te_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)
    k = kc_equilateral(side, m, n)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(k, β, f, μᵣ, εᵣ)
    return te_equilateral_fields(x, y, side, m, n, c_e, c_h)
end

te_equilateral_fields(x::AbstractArray, y::AbstractArray, side, m, n, f, μᵣ, εᵣ) = array_fields(te_equilateral_fields, x, y, side, m, n, f, μᵣ, εᵣ)

"""
    tm_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)

Electric and magnetic field components of an equilateral-triangle TM mode.
Returns `(Ex, Ey, Ez, Hx, Hy, Hz)`.
"""
tm_equilateral_fields(x, y, side, m, n, c_e, c_h) = tm_reflection_fields_from_modal(equilateral_modal_f(x, y, side, m, n, :TM, :S)..., c_e, c_h)

function tm_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)
    k = kc_equilateral(side, m, n)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(k, β, f, μᵣ, εᵣ)
    return tm_equilateral_fields(x, y, side, m, n, c_e, c_h)
end

tm_equilateral_fields(x::AbstractArray, y::AbstractArray, side, m, n, f, μᵣ, εᵣ) = array_fields(tm_equilateral_fields, x, y, side, m, n, f, μᵣ, εᵣ)

"""
    te_right_isosceles_fields(x, y, side, m, n, f, μᵣ, εᵣ)

Electric and magnetic field components of a right-isosceles-triangle TE mode.
Returns `(Ex, Ey, Ez, Hx, Hy, Hz)`.
"""
te_right_isosceles_fields(x, y, side, m, n, c_e, c_h) = te_reflection_fields_from_modal(right_isosceles_modal_f(x, y, side, m, n, :TE)..., c_e, c_h)

function te_right_isosceles_fields(x, y, side, m, n, f, μᵣ, εᵣ)
    k = kc_right_isosceles(side, m, n)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(k, β, f, μᵣ, εᵣ)
    return te_right_isosceles_fields(x, y, side, m, n, c_e, c_h)
end

te_right_isosceles_fields(x::AbstractArray, y::AbstractArray, side, m, n, f, μᵣ, εᵣ)= array_fields(te_right_isosceles_fields, x, y, side, m, n, f, μᵣ, εᵣ)

"""
    tm_right_isosceles_fields(x, y, side, m, n, f, μᵣ, εᵣ)

Electric and magnetic field components of a right-isosceles-triangle TM mode.
Returns `(Ex, Ey, Ez, Hx, Hy, Hz)`.
"""
tm_right_isosceles_fields(x, y, side, m, n, c_e, c_h) = tm_reflection_fields_from_modal(right_isosceles_modal_f(x, y, side, m, n, :TM)..., c_e, c_h)

function tm_right_isosceles_fields(x, y, side, m, n, f, μᵣ, εᵣ)
    k = kc_right_isosceles(side, m, n)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(k, β, f, μᵣ, εᵣ)
    return tm_right_isosceles_fields(x, y, side, m, n, c_e, c_h)
end

tm_right_isosceles_fields(x::AbstractArray, y::AbstractArray, side, m, n, f, μᵣ, εᵣ) = array_fields(tm_right_isosceles_fields, x, y, side, m, n, f, μᵣ, εᵣ)

"""
    te_half_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)

Electric and magnetic field components of a half-equilateral-triangle TE mode.
Returns `(Ex, Ey, Ez, Hx, Hy, Hz)`.
"""
te_half_equilateral_fields(x, y, side, m, n, c_e, c_h) = te_reflection_fields_from_modal(half_equilateral_modal_f(x, y, side, m, n, :TE)..., c_e, c_h)

function te_half_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)
    k = kc_half_equilateral(side, m, n)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = te_coefficients(k, β, f, μᵣ, εᵣ)
    return te_half_equilateral_fields(x, y, side, m, n, c_e, c_h)
end

te_half_equilateral_fields(x::AbstractArray, y::AbstractArray, side, m, n, f, μᵣ, εᵣ) = array_fields(te_half_equilateral_fields, x, y, side, m, n, f, μᵣ, εᵣ)

"""
    tm_half_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)

Electric and magnetic field components of a half-equilateral-triangle TM mode.
Returns `(Ex, Ey, Ez, Hx, Hy, Hz)`.
"""
tm_half_equilateral_fields(x, y, side, m, n, c_e, c_h) = tm_reflection_fields_from_modal(half_equilateral_modal_f(x, y, side, m, n, :TM)..., c_e, c_h)

function tm_half_equilateral_fields(x, y, side, m, n, f, μᵣ, εᵣ)
    k = kc_half_equilateral(side, m, n)
    β = phase_constant(k, f, μᵣ, εᵣ)
    c_e, c_h = tm_coefficients(k, β, f, μᵣ, εᵣ)
    return tm_half_equilateral_fields(x, y, side, m, n, c_e, c_h)
end

tm_half_equilateral_fields(x::AbstractArray, y::AbstractArray, side, m, n, f, μᵣ, εᵣ) = array_fields(tm_half_equilateral_fields, x, y, side, m, n, f, μᵣ, εᵣ)
