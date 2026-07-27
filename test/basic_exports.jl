using Test
using StaticArrays
using LinearAlgebra
using AnalyticEMModes

allfinite_tuple(t) = all(x -> isfinite(real(x)) && isfinite(imag(x)), t)
test_norm(v) = sqrt(sum(abs2, v))
tuple_isapprox(a, b; kwargs...) = all(isapprox(a[i], b[i]; kwargs...) for i in eachindex(a))
tuple_relerr(a, b) = sqrt(sum(abs2, a .- b)) / max(sqrt(sum(abs2, a)), sqrt(sum(abs2, b)), eps())

function spherical_radial_with_derivatives(radial, l, r, k)
    radial == 1 && return AnalyticEMModes.sph_jm_with_derivatives(l, r, k)
    radial == 2 && return AnalyticEMModes.sph_ym_with_derivatives(l, r, k)
    radial == 3 && return AnalyticEMModes.sph_h1m_with_derivatives(l, r, k)
    return AnalyticEMModes.sph_h2m_with_derivatives(l, r, k)
end

function spherical_modal_power(l, m, radial, kind, r, f, μr, εr; Nθ=50, Nϕ=100)
    k = wavenumber(f, μr, εr)
    basis = AnalyticEMModes.SphericalHarmonics(l, normalisation = :sphericart)
    total = 0.0

    for i in 1:Nθ
        θ = (i - 0.5) * π / Nθ
        for j in 1:Nϕ
            ϕ = (j - 0.5) * 2 * π / Nϕ
            x, y, z = AnalyticEMModes.sph2cart(r, θ, ϕ)
            ylm, ∇ylm = AnalyticEMModes.compute_with_gradients(basis, [SVector(x, y, z)])
            idx = l^2 + l + m + 1
            rs = spherical_radial_with_derivatives(radial, l, r, k)
            fields = kind == :TE ?
                     te_from_mn_sph(r, θ, ϕ, rs, ylm[1, idx], ∇ylm[1, idx], l, k, radial, μr, εr) :
                     tm_from_mn_sph(r, θ, ϕ, rs, ylm[1, idx], ∇ylm[1, idx], l, k, radial, μr, εr)
            _, Eθ, Eϕ, _, Hθ, Hϕ = fields
            total += 0.5 * real(Eθ * conj(Hϕ) - Eϕ * conj(Hθ)) * r^2 * sin(θ) * (π / Nθ) * (2 * π / Nϕ)
        end
    end

    return total
end

function spheroidal_pm_angular_check(mode, k, family, even, radial; oblate=false)
    ξ = oblate ? 0.85 : 1.55
    η = 0.27
    ϕ1, ϕ2 = 0.43, 0.91
    q = family == :+ ? mode.m + 1 : mode.m - 1
    v1 = mn_spheroidal_vector(ξ, η, ϕ1, mode, k; family = family, even = even, radial = radial)
    v2 = mn_spheroidal_vector(ξ, η, ϕ2, mode, k; family = family, even = even, radial = radial)
    t12 = even ? (sin(q * ϕ1), sin(q * ϕ2)) : (cos(q * ϕ1), cos(q * ϕ2))
    t345 = even ? (cos(q * ϕ1), cos(q * ϕ2)) : (sin(q * ϕ1), sin(q * ϕ2))
    return tuple_isapprox(v1[[1, 2, 6]] ./ t12[1], v2[[1, 2, 6]] ./ t12[2]; rtol = 1e-9, atol = 1e-8) &&
           tuple_isapprox(v1[[3, 4, 5]] ./ t345[1], v2[[3, 4, 5]] ./ t345[2]; rtol = 1e-9, atol = 1e-8)
end

function spheroidal_curl_m_fd(ξ, η, ϕ, mode, k; family=:z, even=true, radial=1, oblate=false)
    a = real(mode.c / k)
    hξ, hη, hϕ = oblate ? AnalyticEMModes.scale_factors_oblate(a, ξ, η) : AnalyticEMModes.scale_factors_prolate(a, ξ, η)
    sf(u, v) = oblate ? AnalyticEMModes.scale_factors_oblate(a, u, v) : AnalyticEMModes.scale_factors_prolate(a, u, v)
    mvec(u, v, w) = m_spheroidal_vector(u, v, w, mode, k; family = family, even = even, radial = radial, oblate = oblate)
    h = 1e-6
    curlξ = -(((sf(ξ, η + h)[3] * mvec(ξ, η + h, ϕ)[3] - sf(ξ, η - h)[3] * mvec(ξ, η - h, ϕ)[3]) / (2*h) - (hη * mvec(ξ, η, ϕ + h)[2] - hη * mvec(ξ, η, ϕ - h)[2]) / (2*h)) / (hη*hϕ))
    curlη = -(((hξ * mvec(ξ, η, ϕ + h)[1] - hξ * mvec(ξ, η, ϕ - h)[1]) / (2*h) - (sf(ξ + h, η)[3] * mvec(ξ + h, η, ϕ)[3] - sf(ξ - h, η)[3] * mvec(ξ - h, η, ϕ)[3]) / (2*h)) / (hϕ*hξ))
    curlϕ = -(((sf(ξ + h, η)[2] * mvec(ξ + h, η, ϕ)[2] - sf(ξ - h, η)[2] * mvec(ξ - h, η, ϕ)[2]) / (2*h) - (sf(ξ, η + h)[1] * mvec(ξ, η + h, ϕ)[1] - sf(ξ, η - h)[1] * mvec(ξ, η - h, ϕ)[1]) / (2*h)) / (hξ*hη))
    return curlξ, curlη, curlϕ
end

function fd_radial_derivatives(mode, ξ, radial)
    eval_radial = radial == 1 ? AnalyticEMModes.evaluate_radial1 :
                  radial == 2 ? AnalyticEMModes.evaluate_radial2 :
                  radial == 3 ? AnalyticEMModes.evaluate_radial3 :
                  AnalyticEMModes.evaluate_radial4
    h = 1e-6
    R, dR, d²R = eval_radial(mode, ξ)
    R_p, dR_p, _ = eval_radial(mode, ξ + h)
    R_m, dR_m, _ = eval_radial(mode, ξ - h)
    return R, dR, d²R, (R_p - R_m) / (2*h), (dR_p - dR_m) / (2*h)
end

function cartesian_curl_fd(F, x, y, z; h = 1e-6)
    Fx_p, Fx_m = F(x + h, y, z), F(x - h, y, z)
    Fy_p, Fy_m = F(x, y + h, z), F(x, y - h, z)
    Fz_p, Fz_m = F(x, y, z + h), F(x, y, z - h)
    dFdx = (Fx_p .- Fx_m) ./ (2h)
    dFdy = (Fy_p .- Fy_m) ./ (2h)
    dFdz = (Fz_p .- Fz_m) ./ (2h)
    return SVector(dFdy[3] - dFdz[2], dFdz[1] - dFdx[3], dFdx[2] - dFdy[1])
end

function spheroidal_scalar(mode, ξ, η, ϕ; even = true, radial = 1)
    S, _, _ = AnalyticEMModes.evaluate_angular(mode, η)
    R, _, _ = radial == 1 ? AnalyticEMModes.evaluate_radial1(mode, ξ) :
              radial == 2 ? AnalyticEMModes.evaluate_radial2(mode, ξ) :
              radial == 3 ? AnalyticEMModes.evaluate_radial3(mode, ξ) :
                            AnalyticEMModes.evaluate_radial4(mode, ξ)
    return S * R * (even ? cos(mode.m * ϕ) : sin(mode.m * ϕ))
end

function spheroidal_cartesian_vector(mode, family, a, ξ, η, ϕ; oblate = false)
    x, y, z = oblate ? obl2cart(a, ξ, η, ϕ) : pro2cart(a, ξ, η, ϕ)
    if family == :x
        return SVector(1.0, 0.0, 0.0)
    elseif family == :y
        return SVector(0.0, 1.0, 0.0)
    elseif family == :z
        return SVector(0.0, 0.0, 1.0)
    else
        return SVector(x, y, z)
    end
end

function spheroidal_vector_potential(mode, family, a; even = true, radial = 1, oblate = false)
    return (x, y, z) -> begin
        ξ, η, ϕ = oblate ? cart2obl(a, x, y, z) : cart2pro(a, x, y, z)
        ψe = spheroidal_scalar(mode, ξ, η, ϕ; even = true, radial = radial)
        ψo = spheroidal_scalar(mode, ξ, η, ϕ; even = false, radial = radial)
        if family == :+
            return even ? 0.5 * (ψe * SVector(1.0, 0.0, 0.0) - ψo * SVector(0.0, 1.0, 0.0)) :
                          0.5 * (ψo * SVector(1.0, 0.0, 0.0) + ψe * SVector(0.0, 1.0, 0.0))
        elseif family == :-
            return even ? 0.5 * (ψe * SVector(1.0, 0.0, 0.0) + ψo * SVector(0.0, 1.0, 0.0)) :
                          0.5 * (ψo * SVector(1.0, 0.0, 0.0) - ψe * SVector(0.0, 1.0, 0.0))
        end
        spheroidal_scalar(mode, ξ, η, ϕ; even = even, radial = radial) *
            spheroidal_cartesian_vector(mode, family, a, ξ, η, ϕ; oblate = oblate)
    end
end

function mn_spheroidal_cartesian(mode, family, k, ξ, η, ϕ; even = true, radial = 1, oblate = false)
    a = real(mode.c / k)
    mn = mn_spheroidal_vector(ξ, η, ϕ, mode, k; family = family, even = even, radial = radial, oblate = oblate)
    convert_vector = oblate ? oblate_vector_to_cartesian : prolate_vector_to_cartesian
    M = SVector(convert_vector(a, ξ, η, ϕ, mn[1], mn[2], mn[3])...)
    N = SVector(convert_vector(a, ξ, η, ϕ, mn[4], mn[5], mn[6])...)
    return M, N
end

function numerical_mn_spheroidal_cartesian(mode, family, k, ξ, η, ϕ; even = true, radial = 1, oblate = false)
    a = real(mode.c / k)
    x, y, z = oblate ? obl2cart(a, ξ, η, ϕ) : pro2cart(a, ξ, η, ϕ)
    A = spheroidal_vector_potential(mode, family, a; even = even, radial = radial, oblate = oblate)
    M = cartesian_curl_fd(A, x, y, z; h = 1e-6)
    Mfield = (u, v, w) -> cartesian_curl_fd(A, u, v, w; h = 1e-6)
    N = cartesian_curl_fd(Mfield, x, y, z; h = 2e-5) ./ k
    return M, N
end

function spherical_unit_vectors(θ, ϕ)
    er = SVector(sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))
    eθ = SVector(cos(θ) * cos(ϕ), cos(θ) * sin(ϕ), -sin(θ))
    eϕ = SVector(-sin(ϕ), cos(ϕ), 0.0)
    return er, eθ, eϕ
end

function spherical_vector_potential(l, m, k, radial)
    basis = AnalyticEMModes.SphericalHarmonics(l, normalisation = :sphericart)
    idx = l^2 + l + m + 1
    return (x, y, z) -> begin
        r, _, _ = AnalyticEMModes.cart2sph(x, y, z)
        ylm, _ = AnalyticEMModes.compute_with_gradients(basis, [SVector(x, y, z)])
        R, _, _ = spherical_radial_with_derivatives(radial, l, r, k)
        return R * ylm[1, idx] * SVector(x, y, z) / r
    end
end

function mn_spherical_cartesian(l, m, k, radial, p, μr, εr)
    r, θ, ϕ = AnalyticEMModes.cart2sph(p...)
    basis = AnalyticEMModes.SphericalHarmonics(l, normalisation = :sphericart)
    idx = l^2 + l + m + 1
    ylm, ∇ylm = AnalyticEMModes.compute_with_gradients(basis, [p])
    rs = spherical_radial_with_derivatives(radial, l, r, k)
    mn = mn_vectors_sph(r, θ, ϕ, rs, ylm[1, idx], ∇ylm[1, idx], k, μr, εr)
    er, eθ, eϕ = spherical_unit_vectors(θ, ϕ)

    M = mn[1] * er + mn[2] * eθ + mn[3] * eϕ
    N = mn[4] * er + mn[5] * eθ + mn[6] * eϕ
    return M, N
end

function numerical_mn_spherical_cartesian(l, m, k, radial, p)
    A = spherical_vector_potential(l, m, k, radial)
    M = cartesian_curl_fd(A, p[1], p[2], p[3]; h = 1e-6)
    Mfield = (x, y, z) -> cartesian_curl_fd(A, x, y, z; h = 1e-6)
    N = cartesian_curl_fd(Mfield, p[1], p[2], p[3]; h = 2e-5) ./ k
    return M, N
end

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
    @test allfinite_tuple(te[1, 4])
    @test allfinite_tuple(tm[2, 5])
    @test allfinite_tuple(mn[3, 6])

    l, m = 2, 1
    idx = l^2 + l + m + 1
    basis = AnalyticEMModes.SphericalHarmonics(l, normalisation = :sphericart)
    r = test_norm(pts[1])
    rs = AnalyticEMModes.sph_h2m_with_derivatives(l, r, k)
    y, gy = AnalyticEMModes.compute_with_gradients(basis, [pts[1]])
    _, θ, ϕ = AnalyticEMModes.cart2sph(pts[1]...)
    fte = te_from_mn_sph(r, θ, ϕ, rs, y[1, idx], gy[1, idx], l, k, 4, μr, εr)
    ftm = tm_from_mn_sph(r, θ, ϕ, rs, y[1, idx], gy[1, idx], l, k, 4, μr, εr)
    @test allfinite_tuple(fte)
    @test allfinite_tuple(ftm)

    for radial in 1:4, ltest in 0:4, rtest in (0.021, 0.037, 0.064)
        R, dR, d²R_k = spherical_radial_with_derivatives(radial, ltest, rtest, k)
        h = 1e-6
        Rp, dRp, _ = spherical_radial_with_derivatives(radial, ltest, rtest + h, k)
        Rm, dRm, _ = spherical_radial_with_derivatives(radial, ltest, rtest - h, k)
        @test dR ≈ (Rp - Rm) / (2*h) rtol=1e-7 atol=1e-7
        if ltest == 0
            @test abs(d²R_k) < 1e-12
        else
            @test d²R_k ≈ (dRp - dRm) / (2*h) + k^2 * R rtol=1e-6 atol=1e-6
        end
    end

    fd_point = SVector(0.03, 0.011, 0.027)
    for radial in (1, 4), (ltest, mtest) in ((1, 0), (2, 1), (3, -2))
        analyticM, analyticN = mn_spherical_cartesian(ltest, mtest, k, radial, fd_point, μr, εr)
        numericM, numericN = numerical_mn_spherical_cartesian(ltest, mtest, k, radial, fd_point)
        @test tuple_relerr(analyticM, numericM) < 1e-7
        @test tuple_relerr(analyticN, numericN) < 5e-6
    end

    pts2 = [SVector(0.03, 0.011, 0.027), SVector(-0.018, 0.026, 0.041)]
    lmax2 = 3
    basis2 = AnalyticEMModes.SphericalHarmonics(lmax2, normalisation = :sphericart)
    y2, gy2 = AnalyticEMModes.compute_with_gradients(basis2, pts2)

    for radial in (3, 4)
        te_batch = te_sph_fields_lmax(pts2, lmax2, f, μr, εr, radial)
        tm_batch = tm_sph_fields_lmax(pts2, lmax2, f, μr, εr, radial)
        mn_batch = mn_sph_vectors_lmax(pts2, lmax2, f, μr, εr, radial)
        m_batch = m_sph_vectors_lmax(pts2, lmax2, f, μr, εr, radial)
        n_batch = n_sph_vectors_lmax(pts2, lmax2, f, μr, εr, radial)

        @test size(te_batch) == (length(pts2), (lmax2 + 1)^2)
        @test size(tm_batch) == size(te_batch)
        @test size(mn_batch) == size(te_batch)
        @test size(m_batch) == size(te_batch)
        @test size(n_batch) == size(te_batch)

        for pidx in eachindex(pts2), ltest in 1:lmax2, mtest in -ltest:ltest
            idx2 = ltest^2 + ltest + mtest + 1
            r2, θ2, ϕ2 = AnalyticEMModes.cart2sph(pts2[pidx]...)
            rs2 = spherical_radial_with_derivatives(radial, ltest, r2, k)

            te_scalar = te_sph_fields(r2, θ2, ϕ2, rs2, y2[pidx, idx2], gy2[pidx, idx2], k, μr, εr)
            tm_scalar = tm_sph_fields(r2, θ2, ϕ2, rs2, y2[pidx, idx2], gy2[pidx, idx2], k, μr, εr)
            mn_scalar = mn_vectors_sph(r2, θ2, ϕ2, rs2, y2[pidx, idx2], gy2[pidx, idx2], k, μr, εr)
            norm_factor = m_normalization_sph(ltest, test_norm(pts2[1]), k, radial, μr, εr)

            @test tuple_relerr(te_batch[pidx, idx2], te_scalar) < 1e-12
            @test tuple_relerr(tm_batch[pidx, idx2], tm_scalar) < 1e-12
            @test tuple_relerr(mn_batch[pidx, idx2], mn_scalar .* norm_factor) < 1e-12
            @test tuple_isapprox(m_batch[pidx, idx2], mn_batch[pidx, idx2][1:3]; rtol=1e-12, atol=1e-12)
            @test tuple_isapprox(n_batch[pidx, idx2], mn_batch[pidx, idx2][4:6]; rtol=1e-12, atol=1e-12)
        end
    end

    for radial in (3, 4), ltest in 1:3, mtest in -ltest:ltest
        idx2 = ltest^2 + ltest + mtest + 1
        r2, θ2, ϕ2 = AnalyticEMModes.cart2sph(pts2[1]...)
        rs2 = spherical_radial_with_derivatives(radial, ltest, r2, k)
        te_direct = te_sph_fields(r2, θ2, ϕ2, rs2, y2[1, idx2], gy2[1, idx2], k, μr, εr)
        tm_direct = tm_sph_fields(r2, θ2, ϕ2, rs2, y2[1, idx2], gy2[1, idx2], k, μr, εr)
        te_mn = te_from_mn_sph(r2, θ2, ϕ2, rs2, y2[1, idx2], gy2[1, idx2], ltest, k, radial, μr, εr; normalize=false)
        tm_mn = tm_from_mn_sph(r2, θ2, ϕ2, rs2, y2[1, idx2], gy2[1, idx2], ltest, k, radial, μr, εr; normalize=false)
        @test tuple_relerr(te_direct, te_mn) < 1e-12
        @test tuple_relerr(tm_direct, tm_mn) < 1e-12
    end

    for kind in (:TE, :TM), radial in (3, 4), (ltest, mtest) in ((1, 0), (2, 1), (3, -2))
        power = spherical_modal_power(ltest, mtest, radial, kind, 0.043, f, μr, εr)
        expected_power = radial == 3 ? -1.0 : 1.0
        @test power ≈ expected_power rtol=5e-3 atol=5e-3
    end
end

@testset "Basic spheroidal exports" begin
    k = 2π
    @test kc_spheroidal() == 0.0
    @test spheroidal_families() == (:x, :y, :z, :r, :+, :-)

    mode = SpheroidalB(1, 1, 2.4)
    mn = mn_spheroidal_vector(1.4, 0.2, 0.7, mode, k; family = :z, even = true, radial = 4)
    @test allfinite_tuple(mn)

    points = [SVector(1.3, -0.2, 0.1), SVector(1.7, 0.3, 1.1)]
    basis = ProlateSpheroidalBasis(1, 2, 2.4)
    M = mn_spheroidal_vectors(points, basis, k; family = :r, even = false, radial = 4)
    @test size(M) == (length(points), length(basis.basis))
    @test allfinite_tuple(M[1, 1])

    M2 = mn_spheroidal_vectors_mnmax(points, 1, 2, 2.4, k; oblate = false, family = :x, even = true, radial = 4)
    @test size(M2, 1) == length(points)
    @test allfinite_tuple(M2[2, 1])

    prolate_basis = ProlateSpheroidalBasis(2, 3, 2.4)
    oblate_basis = OblateSpheroidalBasis(2, 3, Complex(2.4))
    expected_indices = [(m, n) for m in 0:2 for n in m:3]
    @test [(b.m, b.n) for b in prolate_basis.basis] == expected_indices
    @test [(b.m, b.n) for b in oblate_basis.basis] == expected_indices
    @test all(b -> b.c == 2.4, prolate_basis.basis)
    @test all(b -> b.c == Complex(2.4), oblate_basis.basis)

    ηs = [-0.35, 0.1, 0.42]
    for basis in (prolate_basis, oblate_basis)
        S, dS, d²S = AnalyticEMModes.compute_angular_derivatives(basis, ηs)
        @test size(S) == (length(ηs), length(basis.basis))
        @test size(dS) == size(S)
        @test size(d²S) == size(S)
        for i in eachindex(ηs), j in eachindex(basis.basis)
            scalar = AnalyticEMModes.evaluate_angular(basis.basis[j], ηs[i])
            @test S[i, j] ≈ scalar[1]
            @test dS[i, j] ≈ scalar[2]
            @test d²S[i, j] ≈ scalar[3]
        end
    end

    angular_points = [(0.1, 0.2), (0.42, 1.1)]
    for basis in (prolate_basis, oblate_basis)
        ψ, ∇ψ = AnalyticEMModes.compute_angular_with_derivatives(basis, angular_points)
        @test size(ψ) == (length(angular_points), length(basis.basis))
        @test size(∇ψ) == size(ψ)
        for i in eachindex(angular_points), j in eachindex(basis.basis)
            ηp, ϕp = angular_points[i]
            scalar = AnalyticEMModes.evaluate(basis.basis[j], ηp, ϕp)
            real_scalar = AnalyticEMModes.evaluate_real(basis.basis[j], ηp, ϕp)
            @test ψ[i, j] ≈ scalar[1]
            @test ∇ψ[i, j][1] ≈ scalar[2]
            @test ∇ψ[i, j][2] ≈ scalar[3]
            @test real_scalar[1] ≈ imag(scalar[1])
            @test real_scalar[2] ≈ imag(scalar[2])
            @test real_scalar[3] ≈ imag(scalar[3])
            @test real_scalar[4] ≈ real(scalar[1])
            @test real_scalar[5] ≈ real(scalar[2])
            @test real_scalar[6] ≈ real(scalar[3])
        end
    end

    ξs = [1.25, 1.55]
    oblate_ξs = [0.45, 0.9]
    radial_batches = (
        AnalyticEMModes.compute_radial1_with_derivatives,
        AnalyticEMModes.compute_radial2_with_derivatives,
        AnalyticEMModes.compute_radial3_with_derivatives,
        AnalyticEMModes.compute_radial4_with_derivatives,
    )
    radial_scalars = (
        AnalyticEMModes.evaluate_radial1,
        AnalyticEMModes.evaluate_radial2,
        AnalyticEMModes.evaluate_radial3,
        AnalyticEMModes.evaluate_radial4,
    )
    for (batch_eval, scalar_eval) in zip(radial_batches, radial_scalars), (basis, rs) in ((prolate_basis, ξs), (oblate_basis, oblate_ξs))
        Rb, dRb, d²Rb = batch_eval(basis, rs)
        @test size(Rb) == (length(rs), length(basis.basis))
        @test size(dRb) == size(Rb)
        @test size(d²Rb) == size(Rb)
        for i in eachindex(rs), j in eachindex(basis.basis)
            scalar = scalar_eval(basis.basis[j], rs[i])
            @test Rb[i, j] ≈ scalar[1]
            @test dRb[i, j] ≈ scalar[2]
            @test d²Rb[i, j] ≈ scalar[3]
        end
    end

    oblate_points = [SVector(0.45, -0.2, 0.1), SVector(0.9, 0.3, 1.1)]
    oblate_batch = mn_spheroidal_vectors(oblate_points, oblate_basis, k; family = :z, even = true, radial = 1)
    @test size(oblate_batch) == (length(oblate_points), length(oblate_basis.basis))
    for i in eachindex(oblate_points), j in eachindex(oblate_basis.basis)
        p = oblate_points[i]
        scalar = mn_spheroidal_vector(p[1], p[2], p[3], oblate_basis.basis[j], k; family = :z, even = true, radial = 1)
        @test tuple_isapprox(oblate_batch[i, j], scalar; rtol = 1e-12, atol = 1e-12)
    end

    Mp = mn_spheroidal_vector(1.4, 0.2, 0.7, mode, k; family = :+, even = true, radial = 4)
    Mm = mn_spheroidal_vector(1.4, 0.2, 0.7, mode, k; family = :-, even = false, radial = 4)
    @test allfinite_tuple(Mp)
    @test allfinite_tuple(Mm)

    for oblate in (false, true), even in (true, false), radial in 1:4, family in (:+, :-)
        test_mode = oblate ? SpheroidalB(2, 3, Complex(2.4)) : SpheroidalB(2, 3, 2.4)
        @test spheroidal_pm_angular_check(test_mode, k, family, even, radial; oblate = oblate)
    end

    ξ, η, ϕ = 1.45, 0.27, 0.84
    S, dS, d²S = AnalyticEMModes.evaluate_angular(mode, η)
    R, dR, d²R = AnalyticEMModes.evaluate_radial1(mode, ξ)
    h = 1e-6
    _, dS_p, _ = AnalyticEMModes.evaluate_angular(mode, η + h)
    _, dS_m, _ = AnalyticEMModes.evaluate_angular(mode, η - h)
    _, dR_p, _ = AnalyticEMModes.evaluate_radial1(mode, ξ + h)
    _, dR_m, _ = AnalyticEMModes.evaluate_radial1(mode, ξ - h)
    @test d²S ≈ (dS_p - dS_m) / (2*h) rtol=1e-8 atol=1e-7
    @test d²R ≈ (dR_p - dR_m) / (2*h) rtol=1e-8 atol=1e-7

    oblate_mode = SpheroidalB(1, 1, Complex(2.4))
    _, _, oblate_d²S = AnalyticEMModes.evaluate_angular(oblate_mode, η)
    _, oblate_dS_p, _ = AnalyticEMModes.evaluate_angular(oblate_mode, η + h)
    _, oblate_dS_m, _ = AnalyticEMModes.evaluate_angular(oblate_mode, η - h)
    @test oblate_d²S ≈ (oblate_dS_p - oblate_dS_m) / (2*h) rtol=1e-8 atol=1e-7

    oblate_radial_mode = SpheroidalB(2, 3, Complex(2.4))
    for radial in 1:4
        ξtest = radial == 1 ? 0.85 : 2.5
        _, oblate_dR, oblate_d²R, fd_dR, fd_d²R = fd_radial_derivatives(oblate_radial_mode, ξtest, radial)
        @test oblate_dR ≈ fd_dR rtol=1e-5 atol=1e-5
        @test oblate_d²R ≈ fd_d²R rtol=1e-5 atol=1e-5
    end

    for family in spheroidal_families()
        S, dS, d²S = AnalyticEMModes.evaluate_angular(mode, η)
        R, dR, d²R = AnalyticEMModes.evaluate_radial1(mode, ξ)
        rawM = AnalyticEMModes.eval_M(family, ξ, η, ϕ, S, dS, R, dR, mode.m, mode.c / k, true, :prolate)
        rawN = AnalyticEMModes.eval_N(family, ξ, η, ϕ, S, dS, d²S, R, dR, d²R, mode.m, mode.c / k, k, true, false)
        publicMN = mn_spheroidal_vector(ξ, η, ϕ, mode, k; family = family, even = true, radial = 1)

        @test tuple_isapprox(publicMN[1:3], rawM; rtol = 1e-12, atol = 1e-12)
        @test tuple_isapprox(publicMN[4:6], rawN; rtol = 1e-12, atol = 1e-12)
    end

    for (partials, ξtest) in ((AnalyticEMModes.prolate_partials, ξ), (AnalyticEMModes.oblate_partials, 0.85))
        p = partials(real(mode.c / k), ξtest, η, ϕ)
        eξ = SVector(p[1], p[2], p[3]) / test_norm(SVector(p[1], p[2], p[3]))
        eη = SVector(p[4], p[5], p[6]) / test_norm(SVector(p[4], p[5], p[6]))
        eϕ = SVector(p[7], p[8], p[9]) / test_norm(SVector(p[7], p[8], p[9]))

        @test dot(cross(eξ, eη), eϕ) ≈ -1.0 rtol=1e-12 atol=1e-12
    end

    prolate_numeric_mode = SpheroidalB(2, 3, 2.4)
    oblate_numeric_mode = SpheroidalB(2, 3, Complex(2.4))
    for (test_mode, ξtest, oblate) in ((prolate_numeric_mode, 1.45, false), (oblate_numeric_mode, 0.85, true)),
        family in spheroidal_families(), even in (true, false)

        analyticM, analyticN = mn_spheroidal_cartesian(test_mode, family, k, ξtest, η, ϕ; even = even, radial = 1, oblate = oblate)
        numericM, numericN = numerical_mn_spheroidal_cartesian(test_mode, family, k, ξtest, η, ϕ; even = even, radial = 1, oblate = oblate)

        @test norm(analyticM - numericM) / max(norm(analyticM), norm(numericM), eps()) < 5e-8
        @test norm(analyticN - numericN) / max(norm(analyticN), norm(numericN), eps()) < 2e-4
    end

    for family in spheroidal_families(), even in (true, false)
        curlM = spheroidal_curl_m_fd(ξ, η, ϕ, mode, k; family = family, even = even, radial = 1)
        N = n_spheroidal_vector(ξ, η, ϕ, mode, k; family = family, even = even, radial = 1)
        @test tuple_isapprox(curlM, k .* N; rtol = 2e-8, atol = 2e-7)
    end

    oblate_ξ = 0.85
    for family in spheroidal_families(), even in (true, false)
        curlM = spheroidal_curl_m_fd(oblate_ξ, η, ϕ, oblate_mode, k; family = family, even = even, radial = 1, oblate = true)
        N = n_spheroidal_vector(oblate_ξ, η, ϕ, oblate_mode, k; family = family, even = even, radial = 1, oblate = true)
        @test tuple_isapprox(curlM, k .* N; rtol = 2e-8, atol = 2e-7)
    end
end
