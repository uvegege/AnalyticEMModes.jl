"""
    phase_te(m, k, x)

Phase of the characteristic_coax_equation_te(m, kc, a, b) function.
"""
phase_te(m, k, x) = atan(besselj_prime(m, k * x), bessely_prime(m, k * x)) - atan(besselj_prime(m, x), bessely_prime(m, x))

"""
    phase_tm(m, k, x)

Phase of the characteristic_coax_equation_tm(m, kc, a, b) function.
"""
phase_tm(m, k, x) = atan(besselj(m, k*x), bessely(m, k*x)) - atan(besselj(m, x), bessely(m, x))

"""
    kc_approx_coax_te(m, n, k, t; tol=1.0e-10, max_iters=20)

Returns one of three possible approximations for the zero of the characteristic function of the TE mode.

When `t = 1`, it returns the Cochran approximation (`z¹ps`). 
For `t = 2`, it returns the value of `kc_cwg(k, m, n, :TE)`, while for `t = 3`, it returns the asymptotic approximation (`asymtotic_coax_zero_te`).
The initial estimated value is refined using Newton's method.

"""
function kc_approx_coax_te(m, n, k, t; tol=1.0e-10, max_iters=20)
    if t == 1
        x = z¹ps(m, n-1, k)
    elseif t == 2
        x = kc_cwg(k, m, n, :TE)
    elseif t == 3
        x = asymtotic_coax_zero_te(m, n, 1, k)
    end
    iters = 0
    if x <= 0 || isinf(x)
        return 0.0
    end
    f = characteristic_coax_equation_te(m, x, 1, k)
    while (abs(f) > tol) && (iters < max_iters)
        Δ = f / characteristic_coax_equation_te_prime(m, x, 1, k)
        x -= Δ
        iters += 1
        if x < 0
            return 0.0
        end
        f = characteristic_coax_equation_te(m, x, 1, k)
    end
    return x
end


"""
    verify_approx_te(m, n, k, idmin, vs)

Verify the approximations returned by `kc_approx_coax_te(m, n, k, idmin, vs)`.
Check the convergence of the previous and subsequent zeros (`n-1` and `n+1`) and of a small variation of the ratio `k` (`k-Δk` and `k+Δk`).
"""
function verify_approx_te(m, n, k, idmin, vs)
    # Check n-1, n+1, k-Δk and k+Δk
    vs[idmin] == 0.0 && return false
    v1 = n == 1 ? 0.0 : kc_approx_coax_te(m, n - 1, k, idmin)
    if v1 == 0 && n != 1
        v1 = Inf
    end
    v2 = vs[idmin]
    v3 = kc_approx_coax_te(m, n + 1, k, idmin)
    first_check = (v1 < v2 < v3) && !isapprox(v1, v2) && !isapprox(v1, v3) && !isapprox(v2, v3)
    first_check &= (v2/v1 < 6.0)
    v1 = kc_approx_coax_te(m, n, k + 0.05, idmin)
    v3 = kc_approx_coax_te(m, n, max(k - 0.05, 1.01), idmin)
    second_check = (v1 < v2 < v3) && !isapprox(v1, v2) && !isapprox(v1, v3)  && !isapprox(v2, v3)
    checks = first_check & second_check
    return checks
end


"""
    unknown_approx_te(m, n, k)

Compute a somehow reliable **initial approximation** for the `nth` zero of the Transverse Electric (TE) characteristic equation.

This function employs a robust set of heuristics and criteria designed to overcome the
numerical difficulty of isolating the correct `nth` zero for Bessel functions, 
particularly in complex coaxial or high-order modes.

# Arguments
- `m`: The azimuthal mode number.
- `n`: The radial mode number (the specific zero being sought).
- `k`: The ratio of the outer and inner radii (outer_rad/inner_rad).

# Returns
- `ρ_mn`: The best estimate for the `nth` zero, which can be used as the kc value or as the 
  starting point for an iterative root-finding algorithm.

# Algorithm and Robustness
The procedure is based on comparing the phase characteristics and convergence properties of 
three known asymptotic approximations:

1.  **Candidate Generation:** Three possible approximations are calculated based on known asymptotic series.
2.  **Phase Filtering:** Inaccurate candidates are discarded using strict phase criteria derived from the characteristic function's behavior (see `phase_te`).
3.  **Local Convergence Verification:** Candidates are refined and verified using `verify_approx_te` to ensure they are close to an actual zero.
4.  **Final Heuristics:** A final selection is made based on convergence status and mode index properties.

# Important Note for Advanced Use
While this function is reliable for mode initialization, if unexpected mode ordering are encountered, users are advised to combine their own root-finding 
methods (e.g., using `Roots.jl`) with the characteristic equations:
`characteristic_coax_equation_te` and `characteristic_coax_equation_tm`.
"""
function unknown_approx_te(m, n, k)

    xs1 = z¹ps(m, n-1, k)
    xs2 = kc_cwg(k, m, n, :TE)
    xs3 = asymtotic_coax_zero_te(m, n, 1, k)

    vs1 = kc_approx_coax_te(m, n, k, 1)
    vs2 = kc_approx_coax_te(m, n, k, 2)
    vs3 = kc_approx_coax_te(m, n, k, 3)

    target = iseven(n) ? -1.0 : 1.0

    ct1 = isapprox(cos(phase_te(m, k, vs1)), target, atol=1e-5) && (vs1 > 0.0)
    ct2 = isapprox(cos(phase_te(m, k, vs2)), target, atol=1e-5) && (vs2 > 0.0)
    ct3 = isapprox(cos(phase_te(m, k, vs3)), target, atol=1e-5) && (vs3 > 0.0)
    
    cti = (ct1, ct2, ct3)
    vs = (vs1, vs2, vs3)
    xs = (xs1, xs2, xs3)

    # New
    if count(!=(0.0), vs) == 1 
        id = vs[1] != 0.0 ? 1 : 
        vs[2] != 0.0 ? 2 : 3
        return vs[id]
    end

    if count(cti) == 1
        id = cti[1] ? 1 : 
            cti[2] ? 2 : 3
        return vs[id]
    end

    cti2 = ntuple(id -> begin
            ver = false
            cti[id] || return false
            ver = verify_approx_te(m, n, k, id, vs)
            return ver
        end, 3)

        if count(cti) == 1
        id = cti2[1] ? 1 : 
            cti2[2] ? 2 : 3
        return vs[id]
    end

    ct = ntuple(id -> begin
        if id == 1
            if n > 5 & !(cti2[1] & !cti2[2] & !cti2[3])
                return false
            else
                return cti2[id]
            end
        else
            return cti2[id]
        end
        end, 3)

    if count(ct) == 1
        id = ct[1] ? 1 :
             ct[2] ? 2 : 3
        return vs[id]
    end

    alltrue = true
    for i in 1:3
        for j in 1:3
            i >= j && continue
            (ct[i] & ct[j]) || continue
            alltrue &= isapprox(vs[i], vs[j])
        end
        alltrue || break
    end

    id = ct[1] ? 1 :
         ct[2] ? 2 : 3
    alltrue && return vs[id]

    xmap = map(xs) do xo
        if xo <= 0 | isinf(xo)
            Δ = Inf
        else
            f = characteristic_coax_equation_te(m, xo, 1, k)
            Δ = f / characteristic_coax_equation_te_prime(m, xo, 1, k)
        end
        abs(Δ)
    end

    xmin = Inf
    idmin = 0
    for i in 1:3
        ct[i] || continue
        if xmap[i] < xmin
            xmin = xmap[i]
            idmin = i
        end
    end

    ratios = Inf
    for i in 1:3
        ct[i] || continue
        i == idmin && continue
        isapprox(vs[idmin], vs[i]) && continue
        r = xmap[i] / xmin
        if r < ratios
            ratios = r
        end
    end

    vi1 = isapprox(vs1, vs2)
    vi2 = isapprox(vs1, vs3)
    vi3 = isapprox(vs2, vs3)
    vi = (vi1, vi2, vi3)
    ij = ((1, 2), (1, 3), (2, 3))

    vitrue = false
    for i in 1:3
        ct[i] || continue
        vitrue |= vi[i]
    end

    if vitrue
        conditions = false
        for (idx, v) in enumerate(vi)
            ct[idx] || continue
            if v
                conditions |= min(xmap[ij[idx][1]], xmap[ij[idx][2]]) < 0.025
                conditions && break
            end
        end
    else
        conditions = false
    end

    if conditions
        for (idx, v) in enumerate(vi)
            v == true && return vs[ij[idx][1]]
        end
    end


    if ct[idmin]
        return vs[idmin]
    else
        oidmin = idmin
        omin = Inf
        for i in 1:3
            ct[i] || continue
            i == idmin && continue
            if vs[i] < omin
                oidmin = i
                omin = vs[i]
            end
        end
        return vs[oidmin]
    end
end

"""
    kc_coax_te(a, b, m, n)

Compute the cutoff wavenumber `kc` for the **TE_{m,n}** mode in a coaxial waveguide
using a set of pre-determined analytic approximations based on mode indices (m, n) 
and the radius ratio (k = a/b).

This function acts as a **heuristic selector** to choose the most accurate known
asymptotic formula (`kc_approx_coax_te`) based on an analysis of the problem domain (m, n, k).

# Arguments
- `a`: The outer radius of the coaxial guide.
- `b`: The inner radius of the coaxial guide.
- `m`: The azimuthal mode number.
- `n`: The radial mode number.

# Returns
- `kc`: The approximated cutoff wavenumber.

# Method and Heuristics
The function applies a sequence of criteria (based on observed convergence properties)
to select between different established asymptotic approximations (`kc_approx_coax_te`):
1.  **Special Cases:** Handles the m=0 case and checks for small radius ratios (k < 1.6).
2.  **General Domain:** Uses approximations tailored for lower-order modes.
3.  **High-Ratio Domains:** Switches to different asymptotic series for high radius ratios k >= 15.

If none of the specific criteria are met for the given mode indices and ratio, the function 
defaults to calling `unknown_approx_te(m, n, k)`.
"""
function kc_coax_te(a, b, m, n)
    maxv = max(a, b)
    minv = min(a,b)
    k = maxv/minv
    m == 0 && return kc_approx_coax_te(m, n+1, k, 1)
	k < 1.6 && return kc_approx_coax_te(m, n, k, 1)
	n == 1 && k <= 60 && return kc_approx_coax_te(m, n, k, 1)
	m <= 6 && n <= 2 && return kc_approx_coax_te(m, n, k, 1) # TODO: Test K > 30
	m == 1 && k < 23.5 && return kc_approx_coax_te(m, n, k, 1)
    #m == 0 && return kc_approx_coax_te(m, n+1, k, 1)
	m > 1 && k >= 15 && return kc_approx_coax_te(m, n, k, 2)
	return unknown_approx_te(m, n, k)
end


"""
    kc_approx_coax_tm(m, n, k, t; tol=1.0e-10, max_iters=20)

Returns one of three possible approximations for the zero of the characteristic function of the TE mode.

When `t = 1`, it returns the Cochran approximation (`zps`). 
For `t = 2`, it returns the value of `kc_cwg(k, m, n, :TM)`, while for t = 3, it returns the asymptotic approximation (`asymtotic_coax_zero_tm`).
The initial estimated value is refined using Newton's method.

"""
function kc_approx_coax_tm(m, n, k, t; tol = 1.0e-10, max_iters = 20)
    if t == 1
        x = zps(m, n, k)
    elseif t == 2
        x = kc_cwg(k, m, n, :TM)
    elseif t == 3
        x = asymtotic_coax_zero_tm(m, n, 1, k)
    end
    iters = 0
    if x < 0
        return 0.0
    end
    f = characteristic_coax_equation_tm(m, x, 1, k)
    while (abs(f) > tol) && (iters < max_iters)
        Δ = f / characteristic_coax_equation_tm_prime(m, x, 1, k)
        x -= Δ
        iters += 1
        if x < 0
            return 0.0
        end
        f = characteristic_coax_equation_tm(m, x, 1, k)
    end
    return x
end


"""
    verify_approx_tm(m, n, k, idmin, vs)

Verify the approximations returned by `kc_approx_coax_tm(m, n, k, idmin, vs)`.
Check the convergence of the previous and subsequent zeros (`n-1` and `n+1`) and of a small variation of the ratio `k` (`k-Δk` and `k+Δk`).
"""
function verify_approx_tm(m, n, k, idmin, vs)
    # Check n-1, n+1, k-Δk and k+Δk
    vs[idmin] == 0.0 && return false
    v1 = n  == 1 ? 0.0 : kc_approx_coax_tm(m, n-1, k, idmin)
    v2 = vs[idmin]
    v3 =  kc_approx_coax_tm(m, n+1, k, idmin)
    first_check = (v1 < v2 < v3) && !isapprox(v1, v2) && !isapprox(v1, v3) && !isapprox(v2, v3)
    v1 = kc_approx_coax_tm(m, n, k+0.05, idmin)
    v3 =  kc_approx_coax_tm(m, n, k-0.05, idmin)
    second_check =  (v1 < v2 < v3) && !isapprox(v1, v2) && !isapprox(v1, v3) && !isapprox(v2, v3)
    checks = first_check & second_check
    return checks
end


"""
    unknown_approx_tm(m, n, k)

Calculates a somehow reliable **initial approximation** for the `nth` zero of the Transverse Magnetic (TM) characteristic equation.

This function employs a robust set of heuristics and criteria designed to overcome the
numerical difficulty of isolating the correct `nth` zero for Bessel functions, 
particularly in complex coaxial or high-order modes.

# Arguments
- `m`: The azimuthal mode number.
- `n`: The radial mode number (the specific zero being sought).
- `k`: The ratio of the outer and inner radii (outer_rad/inner_rad).

# Returns
- `ρ_mn`: The best estimate for the `nth` zero, which can be used as the kc value or as the 
  starting point for an iterative root-finding algorithm.

# Algorithm and Robustness
The procedure is based on comparing the phase characteristics and convergence properties of 
three known asymptotic approximations:

1.  **Candidate Generation:** Three possible approximations are calculated based on known asymptotic series.
2.  **Phase Filtering:** Inaccurate candidates are discarded using strict phase criteria derived from the characteristic function's behavior (see `phase_tm`).
3.  **Local Convergence Verification:** Candidates are refined and verified using `verify_approx_tm` to ensure they are close to an actual zero.
4.  **Final Heuristics:** A final selection is made based on convergence status and mode index properties.

# Important Note for Advanced Use
While this function is highly reliable for mode initialization, if numerical stability issues or
unexpected mode ordering are encountered, users are advised to combine their own root-finding 
methods (e.g., using `Roots.jl`) with the rigorous characteristic equations:
`characteristic_coax_equation_te` and `characteristic_coax_equation_tm`.
"""
function unknown_approx_tm(m, n, k)

    xs1 = zps(m, n, k)
    xs2 = kc_cwg(k, m, n, :TM)
    xs3 = asymtotic_coax_zero_tm(m, n, 1, k)

    vs1 = kc_approx_coax_tm(m, n, k, 1)
    vs2 = kc_approx_coax_tm(m, n, k, 2)
    vs3 = kc_approx_coax_tm(m, n, k, 3)

    target = iseven(n) ? 1.0 : -1.0
    ct1 = isapprox(cos(phase_tm(m, k, vs1)), target, atol=1e-5) && (vs1 > 0.0)
    ct2 = isapprox(cos(phase_tm(m, k, vs2)), target, atol=1e-5) && (vs2 > 0.0)
    ct3 = isapprox(cos(phase_tm(m, k, vs3)), target, atol=1e-5) && (vs3 > 0.0)

    
    cti = (ct1, ct2, ct3)
    vs = (vs1, vs2, vs3)
    xs = (xs1, xs2, xs3)

    if count(cti) == 1
        id = cti[1] ? 1 : 
            cti[2] ? 2 : 3
        return vs[id]
    end

    cti2 = ntuple(id -> begin
            ver = false
            cti[id] || return false
            ver = verify_approx_tm(m, n, k, id, vs)
            return ver
        end, 3)

    if count(cti) == 1
        id = cti2[1] ? 1 : 
            cti2[2] ? 2 : 3
        return vs[id]
    end

    ct = ntuple(id -> begin
        if id == 1
            if n > 5 & !(cti2[1] & !cti2[2] & !cti2[3])
                return false
            else
                return cti2[id]
            end
        else
            return cti2[id]
        end
        end, 3)


    if count(ct) == 1
        id = ct[1] ? 1 :
             ct[2] ? 2 : 3
        return vs[id]
    end

    alltrue = true
    for i in 1:3
        for j in 1:3
            i >= j && continue
            (ct[i] & ct[j]) || continue
            alltrue &= isapprox(vs[i], vs[j])
        end
        alltrue || break
    end

    id = ct[1] ? 1 :
         ct[2] ? 2 : 3
    alltrue && return vs[id]

    xmap = map(xs) do xo
        if xo <= 0 | isinf(xo)
            Δ = Inf
        else
            f = characteristic_coax_equation_tm(m, xo, 1, k)
            Δ = f / characteristic_coax_equation_tm_prime(m, xo, 1, k)
        end
        abs(Δ)
    end

    xmin = Inf
    idmin = 0
    for i in 1:3
        ct[i] || continue
        if xmap[i] < xmin
            xmin = xmap[i]
            idmin = i
        end
    end

    ratios = Inf
    for i in 1:3
        ct[i] || continue
        i == idmin && continue
        isapprox(vs[idmin], vs[i]) && continue
        r = xmap[i] / xmin
        if r < ratios
            ratios = r
        end
    end

    vi1 = isapprox(vs1, vs2)
    vi2 = isapprox(vs1, vs3)
    vi3 = isapprox(vs2, vs3)
    vi = (vi1, vi2, vi3)
    ij = ((1, 2), (1, 3), (2, 3))

    vitrue = false
    for i in 1:3
        ct[i] || continue
        vitrue |= vi[i]
    end

    if vitrue
        conditions = false
        for (idx, v) in enumerate(vi)
            ct[idx] || continue
            if v
                conditions |= min(xmap[ij[idx][1]], xmap[ij[idx][2]]) < 0.025
                conditions && break
            end
        end
    else
        conditions = false
    end

    if conditions
        for (idx, v) in enumerate(vi)
            v == true && return vs[ij[idx][1]]
        end
    end


    if ct[idmin]
        return vs[idmin]
    else
        oidmin = idmin
        omin = Inf
        for i in 1:3
            ct[i] || continue
            i == idmin && continue
            if vs[i] < omin
                oidmin = i
                omin = vs[i]
            end
        end
        return vs[oidmin]
    end
end


"""
    kc_coax_tm(a, b, m, n)

Calculates the cutoff wavenumber `kc` for the **TM_{m,n}** mode in a coaxial waveguide
using a set of pre-determined analytic approximations based on mode indices (m, n) 
and the radius ratio (k = a/b).

This function acts as a **heuristic selector** to choose the most accurate known
asymptotic formula (`kc_approx_coax_tm`) based on an analysis of the problem domain (m, n, k).

# Arguments
- `a`: The outer radius of the coaxial guide.
- `b`: The inner radius of the coaxial guide.
- `m`: The azimuthal mode number.
- `n`: The radial mode number.

# Returns
- `kc`: The approximated cutoff wavenumber.

# Method and Heuristics
The function applies a sequence of criteria (based on observed convergence properties)
to select between different established asymptotic approximations (`kc_approx_coax_tm`):
1.  **Special Cases:** Handles the m=0 case and checks for small radius ratios.
2.  **General Domain:** Uses approximations tailored for lower-order modes.
3.  **High-Ratio Domains:** Switches to different asymptotic series for high radius ratios k >= 15.

If none of the specific criteria are met for the given mode indices and ratio, the function 
defaults to calling `unknown_approx_tm(m, n, k)`.
"""
function kc_coax_tm(a, b, m, n)
    maxv = max(a, b)
    minv = min(a,b)
    k = maxv/minv
	k <= 1.5 && return kc_approx_coax_tm(m, n, k, 1) 
	if m < 3 # small m
        if m < 2
            return kc_approx_coax_tm(m, n, k, 1)
        else
            if k < 30.0 
                return kc_approx_coax_tm(m, n, k, 1)
            else		
                return kc_approx_coax_tm(m, n, k, 2)
            end
        end
    end
	k >= 16.0 && kc_approx_coax_tm(m, n, k, 2)
	m >= 11 && k >= 8.0 && kc_approx_coax_tm(m, n, k, 2)
	(n == 1) && (1.5 < k < 16.0) && return kc_approx_coax_tm(m, n, k, 2) 
	(3 <= m <= 6 && n >= 8.0) && return kc_approx_coax_tm(m, n, k, 3)
	return unknown_approx_tm(m, n, k)
end
