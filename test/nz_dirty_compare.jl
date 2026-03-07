include(joinpath(@__DIR__, "..", "src", "spheroidal_vectors.jl"))
include(joinpath(@__DIR__, "helpers", "spheroidal_vectors_diagnostics.jl"))
#include(joinpath(@__DIR__, "..", "src", "spheroidal.jl"))
using LinearAlgebra
using StaticArrays

fit_alpha_complex(h, g) = sum(h .* conj.(g)) / sum(abs2, g)
shape_similarity(h, g) = abs(sum(h .* conj.(g))) / sqrt(sum(abs2, h) * sum(abs2, g))

function run_nz_dirty_compare(; m=2, n=2, c=2.4, k=1.0, even=true, oblate=false, phi=0.7, eta0=0.3,
    hxi=1e-5, heta=1e-5, hphi=1e-5, order=4, verbose=true)
    type = oblate ? :oblate : :prolate
    d = c / k
    basis = oblate ? SpheroidalB(m, n, complex(0.0, c)) : SpheroidalB(m, n, c)

    Sfun(eta) = begin
        S, _, _ = evaluate_angular(basis, eta)
        S
    end
    dSfun(eta) = begin
        _, dS, _ = evaluate_angular(basis, eta)
        dS
    end
    Rfun(xi) = begin
        R, _, _ = evaluate_radial4(basis, xi)
        R
    end
    dRfun(xi) = begin
        _, dR, _ = evaluate_radial4(basis, xi)
        dR
    end

    xis = oblate ? range(0.05, 2.0, length=100) : range(1.05, 3.0, length=100)

    hξ = ComplexF64[]; gsξ = ComplexF64[]; gzξ = ComplexF64[]; grξ = ComplexF64[]
    hη = ComplexF64[]; gsη = ComplexF64[]; gzη = ComplexF64[]; grη = ComplexF64[]
    hφ = ComplexF64[]; gsφ = ComplexF64[]; gzφ = ComplexF64[]; grφ = ComplexF64[]

    for xi in xis
        hNξ, hNη, hNφ = Nᶻ(m, c, xi, eta0, phi, k, basis, even, oblate)
        S, dS, _ = evaluate_angular(basis, eta0)
        R, dR, _ = evaluate_radial4(basis, xi)
        gNsξ, gNsη, gNsφ = N_spheroidal(xi, eta0, phi, S, dS, R, dR, m, k, d, even, type)
        gNzξ, gNzη, gNzφ = Nz_from_M_numeric(
            Mz, xi, eta0, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
            hxi=hxi, heta=heta, hphi=hphi, order=order
        )
        gNrξ, gNrη, gNrφ = Nz_from_M_numeric(
            Mr, xi, eta0, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
            hxi=hxi, heta=heta, hphi=hphi, order=order
        )
        push!(hξ, hNξ); push!(gsξ, gNsξ); push!(gzξ, gNzξ); push!(grξ, gNrξ)
        push!(hη, hNη); push!(gsη, gNsη); push!(gzη, gNzη); push!(grη, gNrη)
        push!(hφ, hNφ); push!(gsφ, gNsφ); push!(gzφ, gNzφ); push!(grφ, gNrφ)
    end

    results_s = Dict(
        "ξ" => (alpha=fit_alpha_complex(hξ, gsξ), sim=shape_similarity(hξ, gsξ)),
        "η" => (alpha=fit_alpha_complex(hη, gsη), sim=shape_similarity(hη, gsη)),
        "ϕ" => (alpha=fit_alpha_complex(hφ, gsφ), sim=shape_similarity(hφ, gsφ)),
    )
    results_z = Dict(
        "ξ" => (alpha=fit_alpha_complex(hξ, gzξ), sim=shape_similarity(hξ, gzξ)),
        "η" => (alpha=fit_alpha_complex(hη, gzη), sim=shape_similarity(hη, gzη)),
        "ϕ" => (alpha=fit_alpha_complex(hφ, gzφ), sim=shape_similarity(hφ, gzφ)),
    )
    results_r = Dict(
        "ξ" => (alpha=fit_alpha_complex(hξ, grξ), sim=shape_similarity(hξ, grξ)),
        "η" => (alpha=fit_alpha_complex(hη, grη), sim=shape_similarity(hη, grη)),
        "ϕ" => (alpha=fit_alpha_complex(hφ, grφ), sim=shape_similarity(hφ, grφ)),
    )

    if verbose
        println("Comparing Nᶻ hardcoded vs N_spheroidal")
        for name in ("ξ", "η", "ϕ")
            r = results_s[name]
            println("Comp $name: alpha=", r.alpha, "  sim=", r.sim)
        end
        println("")
        println("Comparing Nᶻ hardcoded vs curl(Mz)/k")
        for name in ("ξ", "η", "ϕ")
            r = results_z[name]
            println("Comp $name: alpha=", r.alpha, "  sim=", r.sim)
        end
        println("")
        println("Comparing Nᶻ hardcoded vs curl(Mr)/k")
        for name in ("ξ", "η", "ϕ")
            r = results_r[name]
            println("Comp $name: alpha=", r.alpha, "  sim=", r.sim)
        end
    end

    return (from_Nsph=results_s, from_Mz=results_z, from_Mr=results_r)
end

function run_nz_h_sweep(; hs=(1e-4, 3e-5, 1e-5, 3e-6), kwargs...)
    println("h sweep (order=4): η component similarity")
    for h in hs
        res = run_nz_dirty_compare(; hxi=h, heta=h, hphi=h, order=4, verbose=false, kwargs...)
        sim_ns = res.from_Nsph["η"].sim
        sim_mz = res.from_Mz["η"].sim
        sim_mr = res.from_Mr["η"].sim
        println("h=", h, "  simη(Nsph)=", sim_ns, "  simη(curlMz)=", sim_mz, "  simη(curlMr)=", sim_mr)
    end
end

@inline function _partials_oracle(a, ξ, η, ϕ; oblate::Bool)
    if oblate
        ex_x, ex_y, ex_z, eη_x, eη_y, eη_z, eϕ_x, eϕ_y, eϕ_z = oblate_partials(a, ξ, η, ϕ)
    else
        ex_x, ex_y, ex_z, eη_x, eη_y, eη_z, eϕ_x, eϕ_y, eϕ_z = prolate_partials(a, ξ, η, ϕ)
    end
    rξ = SVector(ex_x, ex_y, ex_z)
    rη = SVector(eη_x, eη_y, eη_z)
    rϕ = SVector(eϕ_x, eϕ_y, eϕ_z)
    hξ = norm(rξ); hη = norm(rη); hϕ = norm(rϕ)
    eξ = rξ / hξ; eη = rη / hη; eϕ = rϕ / hϕ
    return (hξ, hη, hϕ, eξ, eη, eϕ)
end

@inline function _curl_orth_oracle(hξ, hη, hϕ,
    dξ_hηAη, dη_hξAξ, dξ_hϕAϕ, dη_hϕAϕ,
    dφ_hηAη, dφ_hξAξ)
    cξ = (dη_hϕAϕ - dφ_hηAη) / (hη * hϕ)
    cη = (dφ_hξAξ - dξ_hϕAϕ) / (hϕ * hξ)
    cϕ = (dξ_hηAη - dη_hξAξ) / (hξ * hη)
    return cξ, cη, cϕ
end

@inline function _fd2_oracle(f, x; dx, xmin=-Inf, xmax=Inf)
    xp = min(x + dx, xmax)
    xm = max(x - dx, xmin)
    if xp == xm
        return (f(xp) - f(x)) / max(dx, eps(x))
    end
    return (f(xp) - f(xm)) / (xp - xm)
end

function Nz_operator_reference(m, c, ξ, η, ϕ, k, mode;
    oblate::Bool=false, even::Bool=true, radial_eval=evaluate_radial4, angular_eval=evaluate_angular, rtol_fd=1e-6)
    a = c / k
    zhat = SVector(0.0, 0.0, 1.0)
    sm, cm = sincos(m * ϕ)
    angϕ = even ? cm : sm
    dangϕ = even ? (-m * sm) : (m * cm)
    ξmin = oblate ? 0.0 : 1.0
    ξmax = Inf
    ηmin, ηmax = -1.0, 1.0
    dξ = max(rtol_fd * max(abs(ξ), 1.0), 1e-8)
    dη = max(rtol_fd * max(abs(η), 1.0), 1e-8)
    dϕ = max(rtol_fd, 1e-6)

    function local_A(ξv, ηv, ϕv, angv, dangv)
        hξ, hη, hϕ, eξ, eη, eϕ = _partials_oracle(a, ξv, ηv, ϕv; oblate=oblate)
        R, _, _ = radial_eval(mode, ξv)
        S, _, _ = angular_eval(mode, ηv)
        ψ = R * S * angv
        ψφ = R * S * dangv
        zξ = dot(zhat, eξ); zη = dot(zhat, eη); zϕ = dot(zhat, eϕ)
        Aξ = zξ * ψ; Aη = zη * ψ; Aϕ = zϕ * ψ
        dφ_hξAξ = hξ * zξ * ψφ
        dφ_hηAη = hη * zη * ψφ
        return (hξ, hη, hϕ, eξ, eη, eϕ, Aξ, Aη, Aϕ, dφ_hξAξ, dφ_hηAη)
    end

    function M_local(ξ0, η0, ϕ0, angv, dangv)
        hξ, hη, hϕ, eξ, eη, eϕ, _, _, _, dφ_hξAξ, dφ_hηAη = local_A(ξ0, η0, ϕ0, angv, dangv)
        f_hϕAϕ_ξ(ξv) = (local_A(ξv, η0, ϕ0, angv, dangv)[3] * local_A(ξv, η0, ϕ0, angv, dangv)[9])
        f_hϕAϕ_η(ηv) = (local_A(ξ0, ηv, ϕ0, angv, dangv)[3] * local_A(ξ0, ηv, ϕ0, angv, dangv)[9])
        f_hηAη(ξv) = (local_A(ξv, η0, ϕ0, angv, dangv)[2] * local_A(ξv, η0, ϕ0, angv, dangv)[8])
        f_hξAξ(ηv) = (local_A(ξ0, ηv, ϕ0, angv, dangv)[1] * local_A(ξ0, ηv, ϕ0, angv, dangv)[7])
        dξ_hϕAϕ = _fd2_oracle(f_hϕAϕ_ξ, ξ0; dx=dξ, xmin=ξmin, xmax=ξmax)
        dη_hϕAϕ = _fd2_oracle(f_hϕAϕ_η, η0; dx=dη, xmin=ηmin, xmax=ηmax)
        dξ_hηAη = _fd2_oracle(f_hηAη, ξ0; dx=dξ, xmin=ξmin, xmax=ξmax)
        dη_hξAξ = _fd2_oracle(f_hξAξ, η0; dx=dη, xmin=ηmin, xmax=ηmax)
        Mξ, Mη, Mϕ = _curl_orth_oracle(hξ, hη, hϕ, dξ_hηAη, dη_hξAξ, dξ_hϕAϕ, dη_hϕAϕ, dφ_hηAη, dφ_hξAξ)
        return (hξ, hη, hϕ, eξ, eη, eϕ, Mξ, Mη, Mϕ)
    end

    hξ, hη, hϕ, eξ, eη, eϕ, _, _, _ = M_local(ξ, η, ϕ, angϕ, dangϕ)
    smp, cmp = sincos(m * (ϕ + dϕ))
    smm, cmm = sincos(m * (ϕ - dϕ))
    angp = even ? cmp : smp
    dangp = even ? (-m * smp) : (m * cmp)
    angm = even ? cmm : smm
    dangm = even ? (-m * smm) : (m * cmm)
    hξp, hηp, _, _, _, _, Mξp, Mηp, _ = M_local(ξ, η, ϕ + dϕ, angp, dangp)
    hξm, hηm, _, _, _, _, Mξm, Mηm, _ = M_local(ξ, η, ϕ - dϕ, angm, dangm)
    dφ_hηMη = (hηp * Mηp - hηm * Mηm) / (2dϕ)
    dφ_hξMξ = (hξp * Mξp - hξm * Mξm) / (2dϕ)

    f_hϕMϕ_ξ(ξv) = begin
        _, _, hϕv, _, _, _, _, _, Mϕv = M_local(ξv, η, ϕ, angϕ, dangϕ)
        hϕv * Mϕv
    end
    f_hϕMϕ_η(ηv) = begin
        _, _, hϕv, _, _, _, _, _, Mϕv = M_local(ξ, ηv, ϕ, angϕ, dangϕ)
        hϕv * Mϕv
    end
    f_hηMη(ξv) = begin
        _, hηv, _, _, _, _, _, Mηv, _ = M_local(ξv, η, ϕ, angϕ, dangϕ)
        hηv * Mηv
    end
    f_hξMξ(ηv) = begin
        hξv, _, _, _, _, _, Mξv, _, _ = M_local(ξ, ηv, ϕ, angϕ, dangϕ)
        hξv * Mξv
    end
    dξ_hϕMϕ = _fd2_oracle(f_hϕMϕ_ξ, ξ; dx=dξ, xmin=ξmin, xmax=ξmax)
    dη_hϕMϕ = _fd2_oracle(f_hϕMϕ_η, η; dx=dη, xmin=ηmin, xmax=ηmax)
    dξ_hηMη = _fd2_oracle(f_hηMη, ξ; dx=dξ, xmin=ξmin, xmax=ξmax)
    dη_hξMξ = _fd2_oracle(f_hξMξ, η; dx=dη, xmin=ηmin, xmax=ηmax)

    Nξ, Nη, Nϕ = _curl_orth_oracle(hξ, hη, hϕ, dξ_hηMη, dη_hξMξ, dξ_hϕMϕ, dη_hϕMϕ, dφ_hηMη, dφ_hξMξ)
    Nξ /= k; Nη /= k; Nϕ /= k
    Ncart = Nξ * eξ + Nη * eη + Nϕ * eϕ
    return Ncart[3]
end

function run_nz_scalar_oracle_compare(; m=2, n=2, c=2.4, k=1.0, even=true, oblate=false, phi=0.7, eta0=0.3)
    type = oblate ? :oblate : :prolate
    d = c / k
    a = d
    basis = oblate ? SpheroidalB(m, n, complex(0.0, c)) : SpheroidalB(m, n, c)
    xis = oblate ? range(0.05, 2.0, length=60) : range(1.05, 3.0, length=60)
    hard = ComplexF64[]; nsph = ComplexF64[]; nmz = ComplexF64[]; oracle = ComplexF64[]
    Sfun(eta) = evaluate_angular(basis, eta)[1]
    dSfun(eta) = evaluate_angular(basis, eta)[2]
    Rfun(xi) = evaluate_radial4(basis, xi)[1]
    dRfun(xi) = evaluate_radial4(basis, xi)[2]
    for xi in xis
        Nξh, Nηh, Nϕh = Nᶻ(m, c, xi, eta0, phi, k, basis, even, oblate)
        S, dS, _ = evaluate_angular(basis, eta0)
        R, dR, _ = evaluate_radial4(basis, xi)
        Nξs, Nηs, Nϕs = N_spheroidal(xi, eta0, phi, S, dS, R, dR, m, k, d, even, type)
        Nξm, Nηm, Nϕm = Nz_from_M_numeric(Mz, xi, eta0, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type; order=4)
        _, _, _, eξ, eη, eϕ = _partials_oracle(a, xi, eta0, phi; oblate=oblate)
        push!(hard, Nξh * eξ[3] + Nηh * eη[3] + Nϕh * eϕ[3])
        push!(nsph, Nξs * eξ[3] + Nηs * eη[3] + Nϕs * eϕ[3])
        push!(nmz, Nξm * eξ[3] + Nηm * eη[3] + Nϕm * eϕ[3])
        push!(oracle, Nz_operator_reference(m, c, xi, eta0, phi, k, basis; oblate=oblate, even=even))
    end
    println("Scalar Nz (= N·ẑ) vs oracle:")
    println("hardcoded: alpha=", fit_alpha_complex(hard, oracle), "  sim=", shape_similarity(hard, oracle))
    println("N_spheroidal: alpha=", fit_alpha_complex(nsph, oracle), "  sim=", shape_similarity(nsph, oracle))
    println("curl(Mz)/k: alpha=", fit_alpha_complex(nmz, oracle), "  sim=", shape_similarity(nmz, oracle))
end

"""
Vector oracle for N^z:
1) build A = zhat * ψ in local spheroidal basis,
2) M = curl(A),
3) N = (1/k) curl(M),
4) return (Nξ,Nη,Nϕ) in local basis.
"""
function Nvec_operator_reference(m, c, ξ, η, ϕ, k, mode;
    oblate::Bool=false, even::Bool=true, radial_eval=evaluate_radial4, angular_eval=evaluate_angular,
    rtol_fd=1e-6)
    a = c / k
    zhat = SVector(0.0, 0.0, 1.0)
    sm, cm = sincos(m * ϕ)
    angϕ = even ? cm : sm
    dangϕ = even ? (-m * sm) : (m * cm)
    ξmin = oblate ? 0.0 : 1.0
    ξmax = Inf
    ηmin, ηmax = -1.0, 1.0
    dξ = max(rtol_fd * max(abs(ξ), 1.0), 1e-8)
    dη = max(rtol_fd * max(abs(η), 1.0), 1e-8)
    dϕ = max(rtol_fd, 1e-6)

    function local_A(ξv, ηv, ϕv, angv, dangv)
        hξ, hη, hϕ, eξ, eη, eϕ = _partials_oracle(a, ξv, ηv, ϕv; oblate=oblate)
        R, _, _ = radial_eval(mode, ξv)
        S, _, _ = angular_eval(mode, ηv)
        ψ = R * S * angv
        ψφ = R * S * dangv
        zξ = dot(zhat, eξ); zη = dot(zhat, eη); zϕ = dot(zhat, eϕ)
        Aξ = zξ * ψ; Aη = zη * ψ; Aϕ = zϕ * ψ
        dφ_hξAξ = hξ * zξ * ψφ
        dφ_hηAη = hη * zη * ψφ
        return (hξ, hη, hϕ, Aξ, Aη, Aϕ, dφ_hξAξ, dφ_hηAη)
    end

    function M_local(ξ0, η0, ϕ0, angv, dangv)
        hξ, hη, hϕ, _, _, _, dφ_hξAξ, dφ_hηAη = local_A(ξ0, η0, ϕ0, angv, dangv)
        f_hϕAϕ_ξ(ξv) = begin
            _, _, hϕv, _, _, Aϕv, _, _ = local_A(ξv, η0, ϕ0, angv, dangv)
            hϕv * Aϕv
        end
        f_hϕAϕ_η(ηv) = begin
            _, _, hϕv, _, _, Aϕv, _, _ = local_A(ξ0, ηv, ϕ0, angv, dangv)
            hϕv * Aϕv
        end
        f_hηAη(ξv) = begin
            _, hηv, _, _, Aηv, _, _, _ = local_A(ξv, η0, ϕ0, angv, dangv)
            hηv * Aηv
        end
        f_hξAξ(ηv) = begin
            hξv, _, _, Aξv, _, _, _, _ = local_A(ξ0, ηv, ϕ0, angv, dangv)
            hξv * Aξv
        end
        dξ_hϕAϕ = _fd2_oracle(f_hϕAϕ_ξ, ξ0; dx=dξ, xmin=ξmin, xmax=ξmax)
        dη_hϕAϕ = _fd2_oracle(f_hϕAϕ_η, η0; dx=dη, xmin=ηmin, xmax=ηmax)
        dξ_hηAη = _fd2_oracle(f_hηAη, ξ0; dx=dξ, xmin=ξmin, xmax=ξmax)
        dη_hξAξ = _fd2_oracle(f_hξAξ, η0; dx=dη, xmin=ηmin, xmax=ηmax)
        Mξ, Mη, Mϕ = _curl_orth_oracle(hξ, hη, hϕ, dξ_hηAη, dη_hξAξ, dξ_hϕAϕ, dη_hϕAϕ, dφ_hηAη, dφ_hξAξ)
        return (hξ, hη, hϕ, Mξ, Mη, Mϕ)
    end

    hξ, hη, hϕ, _, _, _ = M_local(ξ, η, ϕ, angϕ, dangϕ)
    smp, cmp = sincos(m * (ϕ + dϕ))
    smm, cmm = sincos(m * (ϕ - dϕ))
    angp = even ? cmp : smp
    dangp = even ? (-m * smp) : (m * cmp)
    angm = even ? cmm : smm
    dangm = even ? (-m * smm) : (m * cmm)
    hξp, hηp, _, Mξp, Mηp, _ = M_local(ξ, η, ϕ + dϕ, angp, dangp)
    hξm, hηm, _, Mξm, Mηm, _ = M_local(ξ, η, ϕ - dϕ, angm, dangm)
    dφ_hηMη = (hηp * Mηp - hηm * Mηm) / (2dϕ)
    dφ_hξMξ = (hξp * Mξp - hξm * Mξm) / (2dϕ)

    f_hϕMϕ_ξ(ξv) = begin
        _, _, hϕv, _, _, Mϕv = M_local(ξv, η, ϕ, angϕ, dangϕ)
        hϕv * Mϕv
    end
    f_hϕMϕ_η(ηv) = begin
        _, _, hϕv, _, _, Mϕv = M_local(ξ, ηv, ϕ, angϕ, dangϕ)
        hϕv * Mϕv
    end
    f_hηMη(ξv) = begin
        _, hηv, _, _, Mηv, _ = M_local(ξv, η, ϕ, angϕ, dangϕ)
        hηv * Mηv
    end
    f_hξMξ(ηv) = begin
        hξv, _, _, Mξv, _, _ = M_local(ξ, ηv, ϕ, angϕ, dangϕ)
        hξv * Mξv
    end
    dξ_hϕMϕ = _fd2_oracle(f_hϕMϕ_ξ, ξ; dx=dξ, xmin=ξmin, xmax=ξmax)
    dη_hϕMϕ = _fd2_oracle(f_hϕMϕ_η, η; dx=dη, xmin=ηmin, xmax=ηmax)
    dξ_hηMη = _fd2_oracle(f_hηMη, ξ; dx=dξ, xmin=ξmin, xmax=ξmax)
    dη_hξMξ = _fd2_oracle(f_hξMξ, η; dx=dη, xmin=ηmin, xmax=ηmax)

    Nξ, Nη, Nϕ = _curl_orth_oracle(hξ, hη, hϕ, dξ_hηMη, dη_hξMξ, dξ_hϕMϕ, dη_hϕMϕ, dφ_hηMη, dφ_hξMξ)
    return Nξ / k, Nη / k, Nϕ / k
end

function run_nz_vector_oracle_compare(; m=2, n=2, c=2.4, k=1.0, even=true, oblate=false, phi=0.7, eta0=0.3)
    basis = oblate ? SpheroidalB(m, n, complex(0.0, c)) : SpheroidalB(m, n, c)
    xis = oblate ? range(0.05, 2.0, length=60) : range(1.05, 3.0, length=60)
    hξ = ComplexF64[]; hη = ComplexF64[]; hϕ = ComplexF64[]
    oξ = ComplexF64[]; oη = ComplexF64[]; oϕ = ComplexF64[]

    for xi in xis
        Nξh, Nηh, Nϕh = Nᶻ(m, c, xi, eta0, phi, k, basis, even, oblate)
        Nξo, Nηo, Nϕo = Nvec_operator_reference(m, c, xi, eta0, phi, k, basis; oblate=oblate, even=even)
        push!(hξ, Nξh); push!(hη, Nηh); push!(hϕ, Nϕh)
        push!(oξ, Nξo); push!(oη, Nηo); push!(oϕ, Nϕo)
    end

    println("Vector N^z components vs vector operator oracle:")
    println("ξ: alpha=", fit_alpha_complex(hξ, oξ), "  sim=", shape_similarity(hξ, oξ))
    println("η: alpha=", fit_alpha_complex(hη, oη), "  sim=", shape_similarity(hη, oη))
    println("ϕ: alpha=", fit_alpha_complex(hϕ, oϕ), "  sim=", shape_similarity(hϕ, oϕ))
end

function run_li_nz_consistency(; m=2, n=2, c=2.4, k=1.0, even=true, oblate=false, phi=0.7,
    Nξ=40, Nη=40, h=1e-5, order=4, atol=1e-12)
    type = oblate ? :oblate : :prolate
    d = c / k
    basis = oblate ? SpheroidalB(m, n, complex(0.0, c)) : SpheroidalB(m, n, c)
    ξs = oblate ? range(0.05, 3.0, length=Nξ) : range(1.05, 3.5, length=Nξ)
    ηs = range(-0.95, 0.95, length=Nη)

    Sfun(eta) = evaluate_angular(basis, eta)[1]
    dSfun(eta) = evaluate_angular(basis, eta)[2]
    Rfun(xi) = evaluate_radial4(basis, xi)[1]
    dRfun(xi) = evaluate_radial4(basis, xi)[2]

    relξ = Float64[]; relη = Float64[]; relϕ = Float64[]
    absξ = Float64[]; absη = Float64[]; absϕ = Float64[]
    hξv = ComplexF64[]; hηv = ComplexF64[]; hϕv = ComplexF64[]
    gξv = ComplexF64[]; gηv = ComplexF64[]; gϕv = ComplexF64[]

    for ξ in ξs, η in ηs
        Nξh, Nηh, Nϕh = Nᶻ(m, c, ξ, η, phi, k, basis, even, oblate)
        Nξg, Nηg, Nϕg = Nz_from_M_numeric(
            Mz, ξ, η, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
            hxi=h, heta=h, hphi=h, order=order
        )

        Δξ = Nξh - Nξg
        Δη = Nηh - Nηg
        Δϕ = Nϕh - Nϕg

        push!(hξv, Nξh); push!(hηv, Nηh); push!(hϕv, Nϕh)
        push!(gξv, Nξg); push!(gηv, Nηg); push!(gϕv, Nϕg)

        push!(absξ, abs(Δξ)); push!(absη, abs(Δη)); push!(absϕ, abs(Δϕ))
        push!(relξ, abs(Δξ) / max(abs(Nξh), abs(Nξg), atol))
        push!(relη, abs(Δη) / max(abs(Nηh), abs(Nηg), atol))
        push!(relϕ, abs(Δϕ) / max(abs(Nϕh), abs(Nϕg), atol))
    end

    stats = (
        ξ = (max_rel=maximum(relξ), rms_rel=sqrt(sum(x -> x^2, relξ) / length(relξ)),
             max_abs=maximum(absξ), rms_abs=sqrt(sum(x -> x^2, absξ) / length(absξ))),
        η = (max_rel=maximum(relη), rms_rel=sqrt(sum(x -> x^2, relη) / length(relη)),
             max_abs=maximum(absη), rms_abs=sqrt(sum(x -> x^2, absη) / length(absη))),
        ϕ = (max_rel=maximum(relϕ), rms_rel=sqrt(sum(x -> x^2, relϕ) / length(relϕ)),
             max_abs=maximum(absϕ), rms_abs=sqrt(sum(x -> x^2, absϕ) / length(absϕ))),
    )

    αξ = sum(hξv .* conj.(gξv)) / sum(abs2, gξv)
    αη = sum(hηv .* conj.(gηv)) / sum(abs2, gηv)
    αϕ = sum(hϕv .* conj.(gϕv)) / sum(abs2, gϕv)

    errξ = sqrt(sum(abs2, hξv .- αξ .* gξv) / sum(abs2, hξv))
    errη = sqrt(sum(abs2, hηv .- αη .* gηv) / sum(abs2, hηv))
    errϕ = sqrt(sum(abs2, hϕv .- αϕ .* gϕv) / sum(abs2, hϕv))

    shape = (
        ξ=(alpha=αξ, rel_shape_err=errξ),
        η=(alpha=αη, rel_shape_err=errη),
        ϕ=(alpha=αϕ, rel_shape_err=errϕ),
    )

    println("Li consistency check: curl(Mz)/k vs Nᶻ hardcoded")
    println("ξ: ", stats.ξ)
    println("η: ", stats.η)
    println("ϕ: ", stats.ϕ)
    println("Scale-aware fit (Nhard ≈ alpha * Ncurl):")
    println("ξ: ", shape.ξ)
    println("η: ", shape.η)
    println("ϕ: ", shape.ϕ)
    return (raw=stats, shape=shape)
end

run_li_nz_consistency()

run_nz_scalar_oracle_compare(c=2.4)
run_nz_scalar_oracle_compare(c=3.4)
run_nz_scalar_oracle_compare(c=5.4)
run_nz_scalar_oracle_compare(c=1.4)
