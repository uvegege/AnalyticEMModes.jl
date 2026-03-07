# Diagnostic and comparison helpers moved from src/spheroidal_vectors.jl
# Kept under test/helpers to keep src focused on API/runtime vectors.
function metric_factors(xi, eta, d, type)

    if type == :oblate
        f1 = xi^2 + 1
        f2 = xi^2 + eta^2
    else
        f1 = xi^2 - 1
        f2 = xi^2 - eta^2
    end

    f3 = 1 - eta^2

    hξ = d * sqrt(f2 / f1)
    hη = d * sqrt(f2 / f3)
    hφ = d * sqrt(f1 * f3)

    return hξ, hη, hφ, f1, f2, f3
end

function gradpsi(xi, eta, phi, S, dS, R, dR, m, hξ, hη, hφ)

    ψ  = S * R
    dξ = dR * S
    dη = R * dS
    dφ = im * m * ψ

    gξ = dξ / hξ
    gη = dη / hη
    gφ = dφ / hφ

    return gξ, gη, gφ, ψ
end

function M_components(gξ, gη, gφ, xi, eta, f1, f2, f3)
    # ẑ in spheroidal orthonormal basis:
    # ẑ = (eta*sqrt(f1/f2))êξ + (xi*sqrt(f3/f2))êη
    zξ = eta * sqrt(f1 / f2)
    zη = xi  * sqrt(f3 / f2)

    # M = ∇ψ × ẑ  (ẑ_φ = 0 so cross-product terms with gη and gξ against ẑ_φ vanish)
    Mξ = -gφ * zη
    Mη =  gφ * zξ
    Mφ =  gξ * zη - gη * zξ

    return Mξ, Mη, Mφ
end

# Curl of M in orthogonal curvilinear coordinates (product rule applied).
# dηhφ, dξhφ, dξhη, dηhξ are the partial derivatives of the scale factors:
#   prolate: hφ = d√(f1*f3), hη = d√(f2/f3), hξ = d√(f2/f1)
# ∂_φ of all scale factors is zero (axisymmetric).
function curl(Mξ, Mη, Mφ,
    dξMξ, dηMξ, dφMξ,
    dξMη, dηMη, dφMη,
    dξMφ, dηMφ, dφMφ,
    hξ, hη, hφ,
    dηhφ, dξhφ, dξhη, dηhξ)

    Cξ = (1 / (hη * hφ)) * (dηhφ * Mφ + hφ * dηMφ - hη * dφMη)
    Cη = (1 / (hφ * hξ)) * (hξ * dφMξ - dξhφ * Mφ - hφ * dξMφ)
    Cφ = (1 / (hξ * hη)) * (dξhη * Mη + hη * dξMη - dηhξ * Mξ - hξ * dηMξ)

    return Cξ, Cη, Cφ
end


function N_from_M(Cξ, Cη, Cφ, k)

    Nξ = Cξ / k
    Nη = Cη / k
    Nφ = Cφ / k

    return Nξ, Nη, Nφ
end

function N_spheroidal(xi, eta, phi, S, dS, R, dR, m, k, d, even, type)

    if type == :oblate
        f1 = xi^2 + 1
        f2 = xi^2 + eta^2
        s  = 1
    else
        f1 = xi^2 - 1
        f2 = xi^2 - eta^2
        s  = -1
    end

    f3 = 1 - eta^2

    T = xi * dS * R + s * eta * S * dR

    smϕ, cmϕ = sincos(m * phi)
    if !even
        smϕ, cmϕ = -cmϕ, smϕ
    end

    Nξ =  sqrt(f3) / (k * d * sqrt(f2)) * T * cmϕ
    Nη = -sqrt(f1) / (k * d * sqrt(f2)) * T * cmϕ
    Nφ = (m / (k * d)) * (eta / sqrt(f3) - xi / sqrt(f1)) * S * R * smϕ

    return Nξ, Nη, Nφ
end

"""
Compare hardcoded `Mz` against generated `M = ∇ψ × ẑ` at one point.

Inputs `S`, `dS`, `R`, `dR` are interpreted as real amplitudes of the
separated solution `ψ = S(η)R(ξ)e^{imφ}`. The generated field is projected
to real parity basis:
- `even=true`: use `Re(M)`
- `even=false`: use `Im(M)` (equivalent to a `sin/cos` phase shift)
"""
function compare_Mz_hardcoded_generated(xi, eta, phi, S, dS, R, dR, m, d, even, type; atol=1e-12)
    hξ, hη, hφ, f1, f2, f3 = metric_factors(xi, eta, d, type)

    # 1) Hardcoded closed form
    hξ_m, hη_m, hφ_m = Mz(xi, eta, phi, S, dS, R, dR, m, d, even, type)

    # 2) Generated from complex ψ = S*R*phase*exp(i m φ)
    # phase = 1 (even) or -im (odd) to match hardcoded parity convention
    gξ, gη, gφ, _ = gradpsi(xi, eta, phi, S, dS, R, dR, m, hξ, hη, hφ)
    phase = even ? one(ComplexF64) : -im
    A = phase * cis(m * phi)
    cξ, cη, cφ = M_components(A * gξ, A * gη, A * gφ, xi, eta, f1, f2, f3)

    # 3) Hardcoded formulas are in real parity basis
    gξ_m = real(cξ)
    gη_m = real(cη)
    gφ_m = real(cφ)

    # 4) Component-wise errors
    Δξ = hξ_m - gξ_m
    Δη = hη_m - gη_m
    Δφ = hφ_m - gφ_m

    relξ = abs(Δξ) / max(abs(hξ_m), abs(gξ_m), atol)
    relη = abs(Δη) / max(abs(hη_m), abs(gη_m), atol)
    relφ = abs(Δφ) / max(abs(hφ_m), abs(gφ_m), atol)

    return (
        hardcoded = (ξ=hξ_m, η=hη_m, ϕ=hφ_m),
        generated = (ξ=gξ_m, η=gη_m, ϕ=gφ_m),
        abs_error = (ξ=abs(Δξ), η=abs(Δη), ϕ=abs(Δφ)),
        rel_error = (ξ=relξ, η=relη, ϕ=relφ),
    )
end

function _extract_component(v, idx::Int)
    if v isa Tuple
        return v[idx]
    end
    return v
end

function _M_from_profile(mfun, xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, d, even, type)
    # Accept either scalar callbacks or tuple-returning callbacks such as:
    # evaluate_angular -> (S, dS, d2S), evaluate_radial4 -> (R, dR, d2R)
    Sraw  = Sfun(eta)
    dSraw = dSfun(eta)
    Rraw  = Rfun(xi)
    dRraw = dRfun(xi)

    S  = _extract_component(Sraw, 1)
    dS = _extract_component(dSraw, dSraw isa Tuple ? 2 : 1)
    R  = _extract_component(Rraw, 1)
    dR = _extract_component(dRraw, dRraw isa Tuple ? 2 : 1)

    return mfun(xi, eta, phi, S, dS, R, dR, m, d, even, type)
end

function _Mz_from_profile(xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, d, even, type)
    return _M_from_profile(Mz, xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, d, even, type)
end

"""
Compute Nz numerically as `curl(Mz)/k` using central finite differences.

`Sfun,dSfun,Rfun,dRfun` define the profile used by `Mz`.
"""
function Nz_from_Mz_numeric(xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
    hxi=1e-5, heta=1e-5, hphi=1e-5, order=4)
    return Nz_from_M_numeric(Mz, xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
        hxi=hxi, heta=heta, hphi=hphi, order=order)
end

"""
Compute N numerically as `curl(M)/k` using central finite differences,
where `mfun` is one of `Mz`, `Mr`, `Mx`, `My`.
"""
function Nz_from_M_numeric(mfun, xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
    hxi=1e-5, heta=1e-5, hphi=1e-5, order=4)
    if order != 2 && order != 4
        throw(ArgumentError("order must be 2 or 4"))
    end

    function d1(f, x, h)
        if order == 2
            return (f(x + h) - f(x - h)) / (2h)
        else
            return (f(x - 2h) - 8f(x - h) + 8f(x + h) - f(x + 2h)) / (12h)
        end
    end

    function eval_all(xiv, etav, phiv)
        Mξ, Mη, Mφ = _M_from_profile(mfun, xiv, etav, phiv, Sfun, dSfun, Rfun, dRfun, m, d, even, type)
        hξ, hη, hφ, _, _, _ = metric_factors(xiv, etav, d, type)
        return Mξ, Mη, Mφ, hξ, hη, hφ
    end

    _, _, _, hξ, hη, hφ = eval_all(xi, eta, phi)

    dη_hφMφ = d1(ηv -> begin
        _, _, Mφv, _, _, hφv = eval_all(xi, ηv, phi)
        hφv * Mφv
    end, eta, heta)
    dφ_hηMη = d1(ϕv -> begin
        _, Mηv, _, _, hηv, _ = eval_all(xi, eta, ϕv)
        hηv * Mηv
    end, phi, hphi)
    Cξ = (dη_hφMφ - dφ_hηMη) / (hη * hφ)

    dφ_hξMξ = d1(ϕv -> begin
        Mξv, _, _, hξv, _, _ = eval_all(xi, eta, ϕv)
        hξv * Mξv
    end, phi, hphi)
    dξ_hφMφ = d1(ξv -> begin
        _, _, Mφv, _, _, hφv = eval_all(ξv, eta, phi)
        hφv * Mφv
    end, xi, hxi)
    Cη = (dφ_hξMξ - dξ_hφMφ) / (hφ * hξ)

    dξ_hηMη = d1(ξv -> begin
        _, Mηv, _, _, hηv, _ = eval_all(ξv, eta, phi)
        hηv * Mηv
    end, xi, hxi)
    dη_hξMξ = d1(ηv -> begin
        Mξv, _, _, hξv, _, _ = eval_all(xi, ηv, phi)
        hξv * Mξv
    end, eta, heta)
    Cφ = (dξ_hηMη - dη_hξMξ) / (hξ * hη)

    return Cξ / k, Cη / k, Cφ / k
end

"""
Compare hardcoded `N_spheroidal` against numeric `curl(Mz)/k`.
"""
function compare_Nz_hardcoded_generated_numeric(xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
    atol=1e-12, hxi=1e-5, heta=1e-5, hphi=1e-5, order=4)

    S  = Sfun(eta)
    dS = dSfun(eta)
    R  = Rfun(xi)
    dR = dRfun(xi)

    hNξ, hNη, hNφ = N_spheroidal(xi, eta, phi, S, dS, R, dR, m, k, d, even, type)
    gNξ, gNη, gNφ = Nz_from_Mz_numeric(
        xi, eta, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, type;
        hxi=hxi, heta=heta, hphi=hphi, order=order
    )

    Δξ = hNξ - gNξ
    Δη = hNη - gNη
    Δφ = hNφ - gNφ

    relξ = abs(Δξ) / max(abs(hNξ), abs(gNξ), atol)
    relη = abs(Δη) / max(abs(hNη), abs(gNη), atol)
    relφ = abs(Δφ) / max(abs(hNφ), abs(gNφ), atol)

    return (
        hardcoded = (ξ=hNξ, η=hNη, ϕ=hNφ),
        generated = (ξ=gNξ, η=gNη, ϕ=gNφ),
        abs_error = (ξ=abs(Δξ), η=abs(Δη), ϕ=abs(Δφ)),
        rel_error = (ξ=relξ, η=relη, ϕ=relφ),
    )
end


  function fit_Neta_terms(; m=2, n=2, c=2.4, k=1.0, even=true, oblate=false, eta=0.3, phi=0.7)
      basis = oblate ? SpheroidalB(m, n, complex(0.0, c)) : SpheroidalB(m, n, c)
      d = c / k
      xis = oblate ? range(0.05, 2.0, length=120) : range(1.05, 3.0, length=120)

    Sfun(η) = begin
        S, _, _ = evaluate_angular(basis, η)
        S
    end

    dSfun(η) = begin
        S, _, _ = evaluate_angular(basis, η)
        S
    end

    Rfun(ξ) = begin
        R, _ , _ = evaluate_radial4(basis, ξ)
        R
    end

    dRfun(ξ) = begin
        _, dR , _ = evaluate_radial4(basis, ξ)
        dR
    end

      X = Matrix{ComplexF64}(undef, length(xis), 3)
      y = Vector{ComplexF64}(undef, length(xis))

      smϕ, cmϕ = sincos(m * phi)
      if !even
          smϕ, cmϕ = -cmϕ, smϕ
      end

      for (i, xi0) in enumerate(xis)
          xi = abs(xi0)
          factor1 = oblate ? xi^2 + 1   : xi^2 - 1
          factor2 = oblate ? xi^2 + eta^2 : xi^2 - eta^2
          factor3 = 1 - eta^2

          R, dR, d2R = evaluate_radial4(basis, xi0)
          S, dS, _   = evaluate_angular(basis, eta)

          v1 = factor1 / factor2
          v2 = xi * v1

          dv1 = (-2*xi*factor1)/(factor2^2) + (2*xi)/factor2
          dv2 = ((oblate ? 1 : -1) + 3*xi^2)/factor2 + (-2*factor1*xi^2)/(factor2^2)

          T1 = eta * dS * (dv1 * R + v1 * dR)
          T2 = S * (dv2 * dR + v2 * d2R)
          T3 = m^2 * xi / (factor3 * factor1) * S * R

          P = 4 * sqrt(factor3) / (k * d^2 * sqrt(factor2))

          _, Neta_gen, _ = Nz_from_Mz_numeric(
              xi0, eta, phi, Sfun, dSfun, Rfun, dRfun, m, k, d, even, oblate ? :oblate : :prolate;
              hxi=1e-5, heta=1e-5, hphi=1e-5
          )

          X[i, :] = [T1, T2, T3]
          y[i]    = Neta_gen / (P * cmϕ)
      end

      β = X \ y
      yhat = X * β
      rel = norm(y - yhat) / norm(y)

      println("Coeficientes ajustados [a1,a2,a3] en: Nη ~ P*(a1*T1 + a2*T2 + a3*T3)*cmϕ")
      @printf("a1 = %.10f%+.10fim\n", real(β[1]), imag(β[1]))
      @printf("a2 = %.10f%+.10fim\n", real(β[2]), imag(β[2]))
      @printf("a3 = %.10f%+.10fim\n", real(β[3]), imag(β[3]))
      @printf("error relativo ajuste = %.3e\n", rel)

      return β, rel
  end

function _fd1_4(f, x, h)
    return (f(x - 2h) - 8f(x - h) + 8f(x + h) - f(x + 2h)) / (12h)
end

function _v_derivs_x_current_exact(xi, eta, oblate::Bool)
    ξ = abs(xi)
    η = eta
    ξ2 = ξ^2
    η2 = η^2

    factor1 = oblate ? (ξ2 + 1) : (ξ2 - 1)
    factor2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    factor3 = 1 - η2
    ∂ηfactor2 = oblate ? 2η : -2η

    v1 = factor1^(3/2) / factor2
    v2 = ξ * sqrt(factor1) / factor2
    v3 = factor3^(3/2) / factor2
    v4 = η * sqrt(factor3) / factor2
    v5 = sqrt(factor3)
    v6 = sqrt(factor1)
    v7 = η / sqrt(factor3)
    v8 = ξ / sqrt(factor1)

    # Current formulas in Nˣ (as implemented)
    dv1_cur = ξ * sqrt(factor1) * (3 * factor2 - 2 * factor1) / (factor2^2)
    dv2_cur = (((factor1 + ξ2) / sqrt(factor1)) * factor2 - 2 * ξ2 * sqrt(factor1)) / (factor2^2)
    dv3_cur = -3 * η * sqrt(factor3) / factor2 - factor3^(3/2) * ∂ηfactor2 / (factor2^2)
    dv4_cur = (sqrt(factor3) - η2 / sqrt(factor3)) / factor2 - η * sqrt(factor3) * ∂ηfactor2 / (factor2^2)
    dv5_cur = -η / sqrt(factor3)
    dv6_cur = ξ / sqrt(factor1)
    dv7_cur = 1 / sqrt(factor3) + η2 / factor3^(3/2)
    dv8_cur = oblate ? (1 / factor1^(3/2)) : (-1 / factor1^(3/2))

    # Exact derivatives
    dv1_ex = ξ * sqrt(factor1) * (3 * factor2 - 2 * factor1) / factor2^2
    dv2_ex = (((factor1 + ξ2) / sqrt(factor1)) * factor2 - 2 * ξ2 * sqrt(factor1)) / factor2^2
    dv3_ex = -3 * η * sqrt(factor3) / factor2 - factor3^(3/2) * ∂ηfactor2 / factor2^2
    dv4_ex = (sqrt(factor3) - η2 / sqrt(factor3)) / factor2 - η * sqrt(factor3) * ∂ηfactor2 / factor2^2
    dv5_ex = -η / sqrt(factor3)
    dv6_ex = ξ / sqrt(factor1)
    dv7_ex = 1 / sqrt(factor3) + η2 / factor3^(3/2)
    dv8_ex = oblate ? (1 / factor1^(3/2)) : (-1 / factor1^(3/2))

    return (
        v=(v1=v1, v2=v2, v3=v3, v4=v4, v5=v5, v6=v6, v7=v7, v8=v8),
        current=(dv1=dv1_cur, dv2=dv2_cur, dv3=dv3_cur, dv4=dv4_cur, dv5=dv5_cur, dv6=dv6_cur, dv7=dv7_cur, dv8=dv8_cur),
        exact=(dv1=dv1_ex, dv2=dv2_ex, dv3=dv3_ex, dv4=dv4_ex, dv5=dv5_ex, dv6=dv6_ex, dv7=dv7_ex, dv8=dv8_ex),
    )
end

"""
Check derivative quality of v1..v8 used in Nˣ transcription.

Returns max errors for:
- current vs exact formulas
- exact vs finite-difference derivatives
"""
function run_v_derivative_checks_x(; oblate::Bool=false, nξ=120, nη=120, hξ=1e-6, hη=1e-6)
    ξs = oblate ? range(0.2, 3.0, length=nξ) : range(1.05, 3.0, length=nξ)
    ηs = range(-0.8, 0.8, length=nη)

    keys = (:dv1, :dv2, :dv3, :dv4, :dv5, :dv6, :dv7, :dv8)
    max_cur_ex = Dict(k => 0.0 for k in keys)
    max_ex_fd = Dict(k => 0.0 for k in keys)

    for ξ in ξs, η in ηs
        vals = _v_derivs_x_current_exact(ξ, η, oblate)
        cur, ex = vals.current, vals.exact

        max_cur_ex[:dv1] = max(max_cur_ex[:dv1], abs(cur.dv1 - ex.dv1))
        max_cur_ex[:dv2] = max(max_cur_ex[:dv2], abs(cur.dv2 - ex.dv2))
        max_cur_ex[:dv3] = max(max_cur_ex[:dv3], abs(cur.dv3 - ex.dv3))
        max_cur_ex[:dv4] = max(max_cur_ex[:dv4], abs(cur.dv4 - ex.dv4))
        max_cur_ex[:dv5] = max(max_cur_ex[:dv5], abs(cur.dv5 - ex.dv5))
        max_cur_ex[:dv6] = max(max_cur_ex[:dv6], abs(cur.dv6 - ex.dv6))
        max_cur_ex[:dv7] = max(max_cur_ex[:dv7], abs(cur.dv7 - ex.dv7))
        max_cur_ex[:dv8] = max(max_cur_ex[:dv8], abs(cur.dv8 - ex.dv8))

        fd_dv1 = _fd1_4(ξv -> _v_derivs_x_current_exact(ξv, η, oblate).v.v1, ξ, hξ)
        fd_dv2 = _fd1_4(ξv -> _v_derivs_x_current_exact(ξv, η, oblate).v.v2, ξ, hξ)
        fd_dv6 = _fd1_4(ξv -> _v_derivs_x_current_exact(ξv, η, oblate).v.v6, ξ, hξ)

        fd_dv3 = _fd1_4(ηv -> _v_derivs_x_current_exact(ξ, ηv, oblate).v.v3, η, hη)
        fd_dv4 = _fd1_4(ηv -> _v_derivs_x_current_exact(ξ, ηv, oblate).v.v4, η, hη)
        fd_dv5 = _fd1_4(ηv -> _v_derivs_x_current_exact(ξ, ηv, oblate).v.v5, η, hη)
        fd_dv7 = _fd1_4(ηv -> _v_derivs_x_current_exact(ξ, ηv, oblate).v.v7, η, hη)
        fd_dv8 = _fd1_4(ξv -> _v_derivs_x_current_exact(ξv, η, oblate).v.v8, ξ, hξ)

        max_ex_fd[:dv1] = max(max_ex_fd[:dv1], abs(ex.dv1 - fd_dv1))
        max_ex_fd[:dv2] = max(max_ex_fd[:dv2], abs(ex.dv2 - fd_dv2))
        max_ex_fd[:dv6] = max(max_ex_fd[:dv6], abs(ex.dv6 - fd_dv6))
        max_ex_fd[:dv3] = max(max_ex_fd[:dv3], abs(ex.dv3 - fd_dv3))
        max_ex_fd[:dv4] = max(max_ex_fd[:dv4], abs(ex.dv4 - fd_dv4))
        max_ex_fd[:dv5] = max(max_ex_fd[:dv5], abs(ex.dv5 - fd_dv5))
        max_ex_fd[:dv7] = max(max_ex_fd[:dv7], abs(ex.dv7 - fd_dv7))
        max_ex_fd[:dv8] = max(max_ex_fd[:dv8], abs(ex.dv8 - fd_dv8))
    end

    println("=== ", oblate ? "oblate" : "prolate", " : v-derivative checks for Nˣ ===")
    println("max|current-exact|:")
    for k in keys
        println("  ", k, " = ", max_cur_ex[k])
    end
    println("max|exact-fd|:")
    for k in keys
        println("  ", k, " = ", max_ex_fd[k])
    end

    return (current_vs_exact=max_cur_ex, exact_vs_fd=max_ex_fd)
end

