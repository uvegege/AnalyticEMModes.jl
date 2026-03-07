using LinearAlgebra
using StaticArrays
using Gmsh
using Test
using Bessels
using MathieuF

include(joinpath(@__DIR__, "..", "src", "common.jl"))
include(joinpath(@__DIR__, "..", "src", "rectangularwg.jl"))
include(joinpath(@__DIR__, "..", "src", "circularwg.jl"))
include(joinpath(@__DIR__, "..", "src", "coaxialwg.jl"))
include(joinpath(@__DIR__, "..", "src", "ellipticalwg.jl"))
include(joinpath(@__DIR__, "..", "src", "radialwg.jl"))
include(joinpath(@__DIR__, "..", "src", "wedgewg.jl"))
include(joinpath(@__DIR__, "..", "src", "spherical.jl"))
include(joinpath(@__DIR__, "..", "src", "spheroidal.jl"))
include(joinpath(@__DIR__, "..", "src", "spheroidal_power.jl"))

const _MESH_DIR = joinpath(@__DIR__, "..", "docs", "src", "Assets", "mesh")
const _F = 10.0e9  # test frequency [Hz]

# Sz = (1/2) Re(E₁ H₂* - E₂ H₁*) — valid for any orthonormal transverse pair (ê₁, ê₂, êz)
_sz(E1, H2, E2, H1) = 0.5 * real(E1 * conj(H2) - E2 * conj(H1))
# Sρ = (1/2) Re(Eϕ Hz* - Ez Hϕ*) — radial power flux in cylindrical/radial/wedge settings
_sr(Eϕ, Hz, Ez, Hϕ) = 0.5 * real(Eϕ * conj(Hz) - Ez * conj(Hϕ))

#name = raw"C:\MisProyecto\Upload\AnalyticEMModes.jl\docs\src\Assets\mesh\spheroidal_prolate.msh "
                                                                                   
#coord, conn = mesh_data(name)

function mesh_data_radial(name)
    Gmsh.gmsh.initialize()
    Gmsh.gmsh.open(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()

    dt = Gmsh.gmsh.model.getPhysicalGroups(2)[1]
    t = Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dt...)[1]
    element_tags, nodeTagsPerElement = gmsh.model.mesh.getElementsByType(2, t)

    tri_connectivity_flat = nodeTagsPerElement
    num_nodes_per_element = 3 # Triangles
    num_elements = length(tri_connectivity_flat) ÷ num_nodes_per_element
    connectivity = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    nodeTags, nodeCoords, _ = Gmsh.gmsh.model.mesh.getNodes()
    coordinates = reshape(nodeCoords, (3, length(nodeTags)))
    
    dt = Gmsh.gmsh.model.getPhysicalGroups(2)[2]
    t = reduce(vcat, Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dt...))
    et = gmsh.model.mesh.getElementsByType.(2, t)
    tri_connectivity_flat = reduce(vcat, getindex.(et, 2))
    num_elements = length(tri_connectivity_flat) ÷ num_nodes_per_element
    connectivity2 = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    Gmsh.gmsh.finalize()

    return coordinates, connectivity', connectivity2'
end


function mesh_data(name)
    Gmsh.gmsh.initialize()
    Gmsh.gmsh.open(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()

    dt = Gmsh.gmsh.model.getPhysicalGroups(2)[1]
    t = Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dt...)[1]
    _, nodeTagsPerElement = gmsh.model.mesh.getElementsByType(2, t)

    tri_connectivity_flat = nodeTagsPerElement
    num_nodes_per_element = 3 # Triangles
    num_elements = length(tri_connectivity_flat) ÷ num_nodes_per_element
    connectivity = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    nodeTags, nodeCoords, _ = Gmsh.gmsh.model.mesh.getNodes()

    coordinates = reshape(nodeCoords, (3, length(nodeTags)))
    Gmsh.gmsh.finalize()
    return coordinates, connectivity'
end

# ---- utilities: read vertex indices for element e robustly ----
@inline function _tri_nodes(conn, e)
    # conn can be:
    #  - Matrix{Int}: (3, Ne) or (Ne, 3)
    #  - Vector{NTuple{3,Int}} or Vector{SVector{3,Int}} or Vector{Vector{Int}}
    if conn isa AbstractMatrix
        if size(conn,1) == 3
            return (conn[1,e], conn[2,e], conn[3,e])
        elseif size(conn,2) == 3
            return (conn[e,1], conn[e,2], conn[e,3])
        else
            error("conn matrix must be 3×Ne or Ne×3")
        end
    else
        ce = conn[e]
        return (ce[1], ce[2], ce[3])
    end
end

@inline function _pt(coord, i)
    # coord is dim×N (dim=2 or 3)
    if size(coord,1) == 2
        return SVector(coord[1,i], coord[2,i], 0.0)
    elseif size(coord,1) == 3
        return SVector(coord[1,i], coord[2,i], coord[3,i])
    else
        error("coord must be 2×N or 3×N")
    end
end

@inline function tri_area(p1, p2, p3)
    return 0.5 * norm(cross(p2-p1, p3-p1))
end

"""
    integrate_mesh_tri(coord, conn; f, quad=:tri1)

Integrate scalar function f(x::SVector{3}) over a triangular surface mesh.
Returns ∑_tri ∫_tri f dS.

quad:
- :tri1  -> 1-point (centroid), weight = area
- :tri3  -> 3-point degree-2 (barycentric), weights = area/3
"""
function integrate_mesh_tri(coord, conn; f::F, quad::Symbol=:tri1) where F
    Ne = conn isa AbstractMatrix ? (size(conn,1)==3 ? size(conn,2) : size(conn,1)) : length(conn)
    acc = zero(f(_pt(coord,1)))  # works for Float/Complex
    for e in 1:Ne
        i1,i2,i3 = _tri_nodes(conn, e)
        p1 = _pt(coord, i1); p2 = _pt(coord, i2); p3 = _pt(coord, i3)
        A = tri_area(p1,p2,p3)
        if A == 0
            continue
        end
        if quad === :tri1
            xc = (p1+p2+p3)/3
            acc += f(xc) * A
        elseif quad === :tri3
            # barycentric 3-pt rule
            # points: (2/3,1/6,1/6) permutations, weights = A/3
            x1 = (2p1 + p2 + p3)/4
            x2 = (p1 + 2p2 + p3)/4
            x3 = (p1 + p2 + 2p3)/4
            w = A/3
            acc += w*(f(x1) + f(x2) + f(x3))
        else
            error("quad must be :tri1 or :tri3")
        end
    end
    return acc
end

# ===========================================================================
# Power normalization tests
#
# Each testset checks that ∫ Sz dA = 1 when the fields are multiplied by the
# normalization factor returned by the corresponding *_normalization_* function.
# The integral is computed numerically over the waveguide cross-section mesh.
#
# Sz = (1/2) Re(E₁ H₂* - E₂ H₁*) is the Poynting z-component in any orthonormal
# transverse basis (êr, êθ), (êξ, êη), or (êx, êy). Since the meshes are in
# Cartesian (x, y), the dA = dx dy element is automatically handled by triangle
# areas, and the scalar Sz is coordinate-invariant.
# ===========================================================================

# ---- Rectangular Waveguide ----
# rectangular_wg1.msh: W = 1.0, H = 0.5
@testset "Power normalization: rectangular WG" begin
    W, H = 1.0, 0.5
    coord, conn = mesh_data(joinpath(_MESH_DIR, "rectangular_wg1.msh"))

    for (m, n, mode) in [(1, 0, :TE), (2, 0, :TE), (1, 1, :TE), (1, 1, :TM), (2, 1, :TM)]
        kc = kc_rwg(W, H, m, n)
        β  = real(phase_constant(kc, _F, 1.0, 1.0))
        if mode == :TE
            N₀ = te_normalization_rwg(W, H, m, n, kc, β, _F, 1.0, 1.0)
            ce, ch = te_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                Ex, Ey, _, Hx, Hy, _ = te_rwg_fields(p[1], p[2], W, H, m, n, ce, ch)
                _sz(Ex, Hy, Ey, Hx)
            end
        else
            N₀ = tm_normalization_rwg(W, H, m, n, kc, β, _F, 1.0, 1.0)
            ce, ch = tm_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                Ex, Ey, _, Hx, Hy, _ = tm_rwg_fields(p[1], p[2], W, H, m, n, ce, ch)
                _sz(Ex, Hy, Ey, Hx)
            end
        end
        P = N₀^2 * integrate_mesh_tri(coord, conn; f = f_sz, quad = :tri3)
        @test isapprox(real(P), 1.0; rtol = 0.01)
    end
end

# ---- Circular Waveguide ----
# circular_wg1.msh: R = 1.0
@testset "Power normalization: circular WG" begin
    R = 1.0
    coord, conn = mesh_data(joinpath(_MESH_DIR, "circular_wg1.msh"))

    for (m, n, mode) in [(1, 1, :TE), (0, 1, :TE), (0, 1, :TM), (1, 1, :TM)]
        kc = kc_cwg(R, m, n, mode)
        β  = real(phase_constant(kc, _F, 1.0, 1.0))
        if mode == :TE
            N₀ = te_normalization_cwg(R, m, kc, β, _F, 1.0, 1.0)
            ce, ch = te_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                r, θ = hypot(p[1], p[2]), atan(p[2], p[1])
                Er, Eθ, _, Hr, Hθ, _ = te_cwg_fields(r, θ, m, kc, ce, ch)
                _sz(Er, Hθ, Eθ, Hr)
            end
        else
            N₀ = tm_normalization_cwg(R, m, kc, β, _F, 1.0, 1.0)
            ce, ch = tm_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                r, θ = hypot(p[1], p[2]), atan(p[2], p[1])
                Er, Eθ, _, Hr, Hθ, _ = tm_cwg_fields(r, θ, m, kc, ce, ch)
                _sz(Er, Hθ, Eθ, Hr)
            end
        end
        P = N₀^2 * integrate_mesh_tri(coord, conn; f = f_sz, quad = :tri3)
        @test isapprox(real(P), 1.0; rtol = 0.02)
    end
end

# ---- Coaxial Waveguide ----
# coax_wg2.msh: R1 = 1.0 (outer), R2 = 0.1 (inner)
@testset "Power normalization: coaxial WG" begin
    R1, R2 = 1.0, 0.1
    coord, conn = mesh_data(joinpath(_MESH_DIR, "coax_wg2.msh"))

    for (m, n, mode) in [(1, 1, :TE), (0, 1, :TE), (0, 1, :TM), (1, 1, :TM)]
        kc = kc_coax(R1, R2, m, n, mode)
        β  = real(phase_constant(kc, _F, 1.0, 1.0))
        if mode == :TE
            N₀ = te_normalization_coax(R1, R2, m, kc, β, _F, 1.0, 1.0)
            Am, Bm = coax_boundary_coeff_te(m, R2, kc)
            ce, ch = te_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                r, θ = hypot(p[1], p[2]), atan(p[2], p[1])
                Er, Eθ, _, Hr, Hθ, _ = te_coax_fields(r, θ, Am, Bm, m, kc, ce, ch)
                _sz(Er, Hθ, Eθ, Hr)
            end
        else
            N₀ = tm_normalization_coax(R1, R2, m, kc, β, _F, 1.0, 1.0)
            Am, Bm = coax_boundary_coeff_tm(m, R2, kc)
            ce, ch = tm_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                r, θ = hypot(p[1], p[2]), atan(p[2], p[1])
                Er, Eθ, _, Hr, Hθ, _ = tm_coax_fields(r, θ, Am, Bm, m, kc, ce, ch)
                _sz(Er, Hθ, Eθ, Hr)
            end
        end
        P = N₀^2 * integrate_mesh_tri(coord, conn; f = f_sz, quad = :tri3)
        @test isapprox(real(P), 1.0; rtol = 0.02)
    end
end

# ---- Elliptic Waveguide ----
# elliptic_wg1.msh: a = 1.0, b = 0.7
@testset "Power normalization: elliptic WG" begin
    a, b = 1.0, 0.7
    ρ = sqrt(a^2 - b^2)
    coord, conn = mesh_data(joinpath(_MESH_DIR, "elliptic_wg1.msh"))

    for (m, n, even, mode) in [(1, 1, true, :TE), (0, 1, true, :TM), (1, 1, false, :TE), (3, 1, true, :TM)]
        kc = kc_ewg(a, b, m, n, even, mode)
        β  = real(phase_constant(kc, _F, 1.0, 1.0))
        q  = (kc * ρ)^2 / 4
        c_char = even ? MathieuCharA(m, q) : MathieuCharB(m, q)
        coeff  = even ? mathieu_a_coeff(m, q, c_char, 100) : mathieu_b_coeff(m, q, c_char, 100)
        if mode == :TE
            N₀ = te_normalization_ewg(a, b, m, even, kc, β, _F, 1.0, 1.0)
            ce, ch = te_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                ξ, η = cart2elliptic(p[1], p[2], a, b)
                Eξ, Eη, _, Hξ, Hη, _ = te_ewg_fields(ξ, η, ρ, m, even, coeff, q, ce, ch)
                _sz(Eξ, Hη, Eη, Hξ)
            end
        else
            N₀ = tm_normalization_ewg(a, b, m, even, kc, β, _F, 1.0, 1.0)
            ce, ch = tm_coefficients(kc, β, _F, 1.0, 1.0)
            f_sz = p -> begin
                ξ, η = cart2elliptic(p[1], p[2], a, b)
                Eξ, Eη, _, Hξ, Hη, _ = tm_ewg_fields(ξ, η, ρ, m, even, coeff, q, ce, ch)
                _sz(Eξ, Hη, Eη, Hξ)
            end
        end
        P = N₀^2 * integrate_mesh_tri(coord, conn; f = f_sz, quad = :tri3)
        @test isapprox(real(P), 1.0; rtol = 0.03)
    end
end

# ---- Radial Waveguide ----
# radial_wg1.msh: domain is the radial side surface (constant r), PEC are top/bottom planes.
@testset "Power normalization: radial WG" begin
    coord, conn, _ = mesh_data_radial(joinpath(_MESH_DIR, "radial_wg1.msh"))
    R = maximum(hypot.(coord[1, :], coord[2, :]))
    H = maximum(coord[3, :]) - minimum(coord[3, :])

    m, n = 1, 1
    Amn, Bmn = 1.0, 0.0
    β = phase_constant_radial(H, n)
    kc = kc_radial(β, _F, 1.0, 1.0)

    # TE
    N₀_te = te_normalization_radial(R, H, m, Amn, Bmn, kc, _F, 1.0, 1.0)
    ce_te, ch_te = te_coefficients(kc, β, _F, 1.0, 1.0)
    f_sr_te = p -> begin
        r, ϕ, z = hypot(p[1], p[2]), atan(p[2], p[1]), p[3]
        _, Eϕ, Ez, _, Hϕ, Hz = te_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, ce_te, ch_te)
        _sr(Eϕ, Hz, Ez, Hϕ)
    end
    P_te = N₀_te^2 * integrate_mesh_tri(coord, conn; f = f_sr_te, quad = :tri3)
    @test isapprox(real(P_te), 1.0; rtol = 0.03)

    # TM
    N₀_tm = tm_normalization_radial(R, H, m, Amn, Bmn, kc, _F, 1.0, 1.0)
    ce_tm, ch_tm = tm_coefficients(kc, β, _F, 1.0, 1.0)
    f_sr_tm = p -> begin
        r, ϕ, z = hypot(p[1], p[2]), atan(p[2], p[1]), p[3]
        _, Eϕ, Ez, _, Hϕ, Hz = tm_radial_fields(r, ϕ, z, m, Amn, Bmn, kc, β, ce_tm, ch_tm)
        _sr(Eϕ, Hz, Ez, Hϕ)
    end
    P_tm = N₀_tm^2 * integrate_mesh_tri(coord, conn; f = f_sr_tm, quad = :tri3)
    @test isapprox(real(P_tm), 1.0; rtol = 0.03)
end

# ---- Wedge Waveguide ----
# wedge_wg1.msh: domain is the outer radial face (constant r), PEC are side/top/bottom faces.
@testset "Power normalization: wedge WG" begin
    coord, conn, _ = mesh_data_radial(joinpath(_MESH_DIR, "wedge_wg1.msh"))
    R = maximum(hypot.(coord[1, :], coord[2, :]))
    H = maximum(coord[3, :]) - minimum(coord[3, :])
    ϕ_vals = atan.(coord[2, :], coord[1, :])
    ϕ0 = maximum(ϕ_vals) - minimum(ϕ_vals)

    p_idx, n = 1, 1
    Amn, Bmn = 1.0, 0.0
    β = phase_constant_radial(H, n)
    kc = kc_wedge(β, _F, 1.0, 1.0)

    # TE
    N₀_te = te_normalization_wedge(R, H, ϕ0, p_idx, Amn, Bmn, kc, _F, 1.0, 1.0)
    ce_te, ch_te = te_coefficients(kc, β, _F, 1.0, 1.0)
    f_sr_te = p -> begin
        r, ϕ, z = hypot(p[1], p[2]), atan(p[2], p[1]), p[3]
        _, Eϕ, Ez, _, Hϕ, Hz = te_wedge_fields(r, ϕ, z, ϕ0, p_idx, Amn, Bmn, kc, β, ce_te, ch_te)
        _sr(Eϕ, Hz, Ez, Hϕ)
    end
    P_te = N₀_te^2 * integrate_mesh_tri(coord, conn; f = f_sr_te, quad = :tri3)
    @test isapprox(real(P_te), 1.0; rtol = 0.04)

    # TM
    N₀_tm = tm_normalization_wedge(R, H, ϕ0, p_idx, Amn, Bmn, kc, _F, 1.0, 1.0)
    ce_tm, ch_tm = tm_coefficients(kc, β, _F, 1.0, 1.0)
    f_sr_tm = p -> begin
        r, ϕ, z = hypot(p[1], p[2]), atan(p[2], p[1]), p[3]
        _, Eϕ, Ez, _, Hϕ, Hz = tm_wedge_fields(r, ϕ, z, ϕ0, p_idx, Amn, Bmn, kc, β, ce_tm, ch_tm)
        _sr(Eϕ, Hz, Ez, Hϕ)
    end
    P_tm = N₀_tm^2 * integrate_mesh_tri(coord, conn; f = f_sr_tm, quad = :tri3)
    @test isapprox(real(P_tm), 1.0; rtol = 0.04)
end

# ---- Spherical Region ----
# sphregion.msh: sphere surface.
@testset "Power normalization: spherical region (M/N)" begin
    coord, conn = mesh_data(joinpath(_MESH_DIR, "sphregion.msh"))
    R0 = mean(hypot.(coord[1, :], coord[2, :], coord[3, :]))
    k = wavenumber(_F, 1.0, 1.0)
    radial = 4
    basis = SphericalHarmonics(1, normalisation = :sphericart)

    # Test l=1,m=0
    l, m = 1, 0
    idx_lm = l^2 + l + m + 1

    # TE
    N₀_te = te_normalization_sph(l, R0, k, radial, 1.0, 1.0)
    f_sr_te = p -> begin
        r, θ, ϕ = cart2sph(p[1], p[2], p[3])
        rs = radial == 3 ? sph_h1m_with_derivatives(l, r, k) : sph_h2m_with_derivatives(l, r, k)
        Y, ∇Y = compute_with_gradients(basis, [p])
        _, Eθ, Eϕ, _, Hθ, Hϕ = te_sph_fields(r, θ, ϕ, rs, Y[1, idx_lm], ∇Y[1, idx_lm], k, 1.0, 1.0)
        _sz(Eθ, Hϕ, Eϕ, Hθ) # radial flux in spherical basis
    end
    P_te = N₀_te^2 * integrate_mesh_tri(coord, conn; f = f_sr_te, quad = :tri3)
    @test isapprox(real(P_te), 1.0; rtol = 0.04)

    # TM
    N₀_tm = tm_normalization_sph(l, R0, k, radial, 1.0, 1.0)
    f_sr_tm = p -> begin
        r, θ, ϕ = cart2sph(p[1], p[2], p[3])
        rs = radial == 3 ? sph_h1m_with_derivatives(l, r, k) : sph_h2m_with_derivatives(l, r, k)
        Y, ∇Y = compute_with_gradients(basis, [p])
        _, Eθ, Eϕ, _, Hθ, Hϕ = tm_sph_fields(r, θ, ϕ, rs, Y[1, idx_lm], ∇Y[1, idx_lm], k, 1.0, 1.0)
        _sz(Eθ, Hϕ, Eϕ, Hθ) # radial flux in spherical basis
    end
    P_tm = N₀_tm^2 * integrate_mesh_tri(coord, conn; f = f_sr_tm, quad = :tri3)
    @test isapprox(real(P_tm), 1.0; rtol = 0.04)
end

# ---- Spheroidal Region ----
# spheroidal_prolate.msh / spheroidal_oblate.msh: closed ξ=const surfaces.
@testset "Power normalization: spheroidal region (M/N)" begin
    k = wavenumber(_F, 1.0, 1.0)
    m, n = 1, 1
    family = :z
    radial = 4

    for (meshname, oblate) in [("spheroidal_prolate.msh", false), ("spheroidal_oblate.msh", true)]
        coord, conn = mesh_data(joinpath(_MESH_DIR, meshname))
        ρmax = maximum(hypot.(coord[1, :], coord[2, :]))
        zmax = maximum(abs.(coord[3, :]))
        major, minor = max(ρmax, zmax), min(ρmax, zmax)
        a = sqrt(max(major^2 - minor^2, eps()))
        c = k * a
        mode = oblate ? SpheroidalB(m, n, Complex(c)) : SpheroidalB(m, n, c)

        # Estimate ξ0 from one surface point
        p0 = SVector(coord[1, 1], coord[2, 1], coord[3, 1])
        ξ0 = if oblate
            ξ, _, _ = cart2obl(a, p0[1], p0[2], p0[3]); ξ
        else
            ξ, _, _ = cart2pro(a, p0[1], p0[2], p0[3]); ξ
        end

        # Normalization from semi-analytical kernel approach
        N₀ = spheroidal_mn_normalization_on_xi_surface_1d(mode, k, ξ0;
                                                          family = family, even = true,
                                                          oblate = oblate, radial = radial,
                                                          CE = 1.0 + 0im, CH = 1.0 + 0im,
                                                          E_from = :M, H_from = :N,
                                                          Nη = 257)

        # Independent mesh-based flux integral over spheroidal surface
        f_sξ = p -> begin
            ξ, η, ϕ = if oblate
                cart2obl(a, p[1], p[2], p[3])
            else
                cart2pro(a, p[1], p[2], p[3])
            end
            Mξ, Mη, Mϕ, Nξ, Nη, Nϕ = mn_spheroidal_vector(ξ, η, ϕ, mode, k;
                                                          family = family, even = true,
                                                          oblate = oblate, radial = radial)
            _sz(Mη, Nϕ, Mϕ, Nη)
        end
        P = N₀^2 * integrate_mesh_tri(coord, conn; f = f_sξ, quad = :tri3)
        @show real(P), N₀
        @test isapprox(real(P), 1.0; rtol = 0.08)
    end
end
