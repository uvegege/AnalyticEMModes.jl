using Test
using StaticArrays
using AnalyticEMModes

_allfinite_tuple(t) = all(x -> isfinite(real(x)) && isfinite(imag(x)), t)
_test_norm(v) = sqrt(sum(abs2, v))

@testset "Basic exported API" begin
    f = 10e9
    μr, εr = 1.0, 1.0

    k = wavenumber(f, μr, εr)
    @test k > 0
    @test cutoff_frequency(2π, μr, εr) > 0
    @test mode_wavelength(2π) ≈ 1.0
    @test mode_impedance(:TE, cutoff_frequency(2π, μr, εr), f, μr, εr) isa Number

    @test length(first_n_modes_rwg(6, 1.0, 0.5)) == 6
    @test length(first_n_modes_cwg(6, 1.0)) == 6
    @test length(first_n_modes_coax(6, 1.0, 0.3)) == 6
    @test length(first_n_modes_radial(6, 0.8)) == 6
    @test length(first_n_modes_ewg(6, 1.0, 0.7)) == 6
    @test length(first_n_modes_elliptic_radial(6, 0.8)) == 6
    @test length(first_n_modes_sph(6)) == 6
end

@testset "Elliptic radial normalization" begin
    focal_distance = 0.35
    height = 0.8
    μr, εr = 1.0, 1.0
    f = 2.0 * cutoff_frequency_elliptic_radial(height, 1, μr, εr)
    ξ = 1.1

    te_mode = AnalyticEMModes.EllipticRadialMode(focal_distance, height, 1, 1, true, :TE, f, μr, εr)
    tm_mode = AnalyticEMModes.EllipticRadialMode(focal_distance, height, 1, 1, false, :TM, f, μr, εr)

    Ate = te_normalization_elliptic_radial(te_mode, ξ, f, μr, εr; Nη=512)
    Atm = tm_normalization_elliptic_radial(tm_mode, ξ, f, μr, εr; Nη=512)

    @test isfinite(Ate) && Ate > 0
    @test isfinite(Atm) && Atm > 0
    @test Ate^2 * AnalyticEMModes.elliptic_radial_modal_power_1d(te_mode, ξ, f, μr, εr; Nη=512) ≈ 1.0
    @test Atm^2 * AnalyticEMModes.elliptic_radial_modal_power_1d(tm_mode, ξ, f, μr, εr; Nη=512) ≈ 1.0
end

@testset "Basic spherical exports" begin
    f = 10e9
    μr, εr = 1.0, 1.0
    k = wavenumber(f, μr, εr)
    pts = [SVector(0.03, 0.0, 0.0), SVector(0.0, 0.03, 0.0), SVector(0.0, 0.0, 0.03)]
    lmax = 2

    te = te_sph_fields_lmax(pts, lmax, f, μr, εr, 4)
    tm = tm_sph_fields_lmax(pts, lmax, f, μr, εr, 4)
    mn = mn_sph_vectors_lmax(pts, lmax, f, μr, εr, 4)
    @test size(te) == (length(pts), (lmax + 1)^2)
    @test size(tm) == (length(pts), (lmax + 1)^2)
    @test size(mn) == (length(pts), (lmax + 1)^2)
    @test _allfinite_tuple(te[1, 4])
    @test _allfinite_tuple(tm[2, 5])
    @test _allfinite_tuple(mn[3, 6])

    l, m = 2, 1
    idx = l^2 + l + m + 1
    basis = AnalyticEMModes.SphericalHarmonics(l, normalisation = :sphericart)
    r = _test_norm(pts[1])
    rs = AnalyticEMModes.sph_h2m_with_derivatives(l, r, k)
    y, gy = AnalyticEMModes.compute_with_gradients(basis, [pts[1]])
    _, θ, ϕ = AnalyticEMModes.cart2sph(pts[1]...)
    fte = te_from_mn_sph(r, θ, ϕ, rs, y[1, idx], gy[1, idx], l, k, 4, μr, εr)
    ftm = tm_from_mn_sph(r, θ, ϕ, rs, y[1, idx], gy[1, idx], l, k, 4, μr, εr)
    @test _allfinite_tuple(fte)
    @test _allfinite_tuple(ftm)
end

@testset "Basic spheroidal exports" begin
    k = 2π
    @test kc_spheroidal() == 0.0
    @test spheroidal_families() == (:x, :y, :z, :r)

    mode = SpheroidalB(1, 1, 2.4)
    mn = mn_spheroidal_vector(1.4, 0.2, 0.7, mode, k; family = :z, even = true, radial = 4)
    @test _allfinite_tuple(mn)

    points = [SVector(1.3, -0.2, 0.1), SVector(1.7, 0.3, 1.1)]
    basis = ProlateSpheroidalBasis(1, 2, 2.4)
    M = mn_spheroidal_vectors(points, basis, k; family = :r, even = false, radial = 4)
    @test size(M) == (length(points), length(basis.basis))
    @test _allfinite_tuple(M[1, 1])

    M2 = mn_spheroidal_vectors_mnmax(points, 1, 2, 2.4, k; oblate = false, family = :x, even = true, radial = 4)
    @test size(M2, 1) == length(points)
    @test _allfinite_tuple(M2[2, 1])
end
