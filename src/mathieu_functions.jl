# Computation of Mathieu angular and radial functions with their respective derivatives.
# References: 
# 1. Algorithms for the Computation of all Mathieu Functions of Integer Orders, Fayez A. Alhargan.
# 2. Implementation of angular and radial Mathieu function based on GSL.

#TODO: Implemente a/b Char value or keep using MathieuF?

# Angular Functions 

"""
    ce_kernel(order, coeff, q, z)

Kernel for computing the Even Mathieu function `ce_m` and its derivative.
"""
function ce_kernel(order, coeff, q, z)

    if q == 0
        norm = ifelse(order == 0, sqrt(2), one(z))
        return (cos(order * z) / norm, -sin(order * z) * order)
    end

    f = zero(z)
    f¹ = zero(z)
    norm = zero(real(f))
    if iseven(order)
        for i in eachindex(coeff)
            norm += coeff[i]^2
            f += coeff[i] * cos(2 * (i - 1) * z)
            f¹ -= coeff[i] * sin(2 * (i - 1) * z) * 2 * (i - 1)
        end
    else
        for i in eachindex(coeff)
            norm += coeff[i]^2
            f += coeff[i] * cos((2 * (i - 1) + 1) * z)
            f¹ -= coeff[i] * sin((2 * (i - 1) + 1) * z) * (2 * (i - 1) + 1)
        end
    end
    f = f / norm
    f¹ = f¹ / norm
    return f, f¹
end

"""
    se_kernel(order, coeff, q, z)

Kernel for computing the Odd Mathieu function `se_m` and its derivative.
"""
function se_kernel(order, coeff, q, z)

    order == 0 && return (zero(z), zero(z))
    q == 0 && return (sin(order * z), order * cos(order * z))

    f = zero(z)
    f¹ = zero(z)
    norm = zero(real(f))
    if iseven(order)
        for i in eachindex(coeff)
            norm += coeff[i]^2
            f += coeff[i] * sin(2 * i * z)
            f¹ += coeff[i] * cos((2 * i) * z) * (2 * i)
        end
    else
        for i in eachindex(coeff)
            norm += coeff[i]^2
            f += coeff[i] * sin((2 * (i - 1) + 1) * z)
            f¹ += coeff[i] * cos((2 * (i - 1) + 1) * z) * (2 * (i - 1) + 1)
        end
    end
    f = f / norm
    f¹ = f¹ / norm
    return f, f¹
end


"""
    mathieu_ce(order, q, z)
    mathieu_ce(order, q, z::AbstractArray)

Computes the Even Mathieu function `ce_m` and its derivative.
Returns the function of order m and parameter q evaluated at z (in rad).

Computes the eigenvalue and coefficients only once when `z` is an array. 

"""
function mathieu_ce(order, q, z)
    (order < 0) && (order = -order)
    a = MathieuCharA(order, q)
    coeff = mathieu_a_coeff(order, q, a, 100)
    f, f¹ = ce_kernel(order, coeff, q, z)
    return f, f¹
end

"""
    mathieu_se(order, q, z)
    mathieu_se(order, q, z::AbstractArray)
    
Computes the Odd Mathieu function `se_m` and its derivative.
Returns the function of order m and parameter q evaluated at z (in rad).

Computes the eigenvalue and coefficients only once when `z` is an array. 

"""
function mathieu_se(order, q, z)
    (order < 0) && (order = -order)
    a = MathieuCharB(order, q)
    coeff = mathieu_b_coeff(order, q, a, 100)
    f, f¹ = se_kernel(order, coeff, q, z)
    return f, f¹
end

function mathieu_ce(order, q, z::AbstractArray)
    (order < 0) && (order = -order)
    a = MathieuCharA(order, q)
    coeff = mathieu_a_coeff(order, q, a, 100)
    f = ce_kernel.(order, Ref(coeff), q, z)
    return f
end

function mathieu_se(order, q, z::AbstractArray)
    (order < 0) && (order = -order)
    a = MathieuCharB(order, q)
    coeff = mathieu_b_coeff(order, q, a, 100)
    f = se_kernel.(order, Ref(coeff), q, z)
    return f
end

# Radial Functions
"""
    Mce_kernel(kind, order, coeff, q, z; maxerr=1e-14)

Kernel for computing the Even Modified Mathieu function `mathieu_Mce` and its derivative of `kind` (1 or 2) of order `m`, parameter `q` at `z`.
"""
function Mce_kernel(kind, order, coeff, q, z; maxerr=1e-14)

    amax = 0.0
    fn = 0.0
    fn¹ = zero(fn)
    u1 = sqrt(q) * exp(-z)
    u2 = sqrt(q) * exp(z)
    du1_dz = -u1
    du2_dz = u2

    if iseven(order)
        for k in 0:length(coeff)-1
            amax = max(amax, abs(coeff[k+1]))
            (abs(coeff[k+1]) / amax < maxerr) && break
            j1c = besselj(k, u1)
            j1c_prime = besselj_prime(k, u1)
            z2c = kind == 1 ? besselj(k, u2) : bessely(k, u2)
            z2c_prime = kind == 1 ? besselj_prime(k, u2) : bessely_prime(k, u2)
            fc = (-1.0)^(0.5 * order + k) * coeff[k+1]
            fn += fc * j1c * z2c
            fn¹ += fc * (j1c * z2c_prime * du2_dz + j1c_prime * z2c * du1_dz)
        end
    else
        for k in 0:length(coeff)-1
            amax = max(amax, abs(coeff[k+1]))
            abs(coeff[k+1]) / amax < maxerr && break
            j1c = besselj(k, u1)
            j1pc = besselj(k + 1, u1)
            j1c_prime = besselj_prime(k, u1)
            j1pc_prime = besselj_prime(k + 1, u1)
            fl = kind == 1
            z2c = fl ? besselj(k, u2) : bessely(k, u2)
            z2pc = fl ? besselj(k + 1, u2) : bessely(k + 1, u2)
            z2c_prime = fl ? besselj_prime(k, u2) : bessely_prime(k, u2)
            z2pc_prime = fl ? besselj_prime(k + 1, u2) : bessely_prime(k + 1, u2)

            fc = (-1.0)^(0.5 * (order - 1) + k) * coeff[k+1]
            fn += fc * (j1c * z2pc + j1pc * z2c)
            fn¹ += fc * (j1c * z2pc_prime * du2_dz + j1c_prime * z2pc * du1_dz +
                         +j1pc * z2c_prime * du2_dz + j1pc_prime * z2c * du1_dz)
        end
    end
    fn *= sqrt(π / 2.0) / coeff[1]
    fn¹ *= sqrt(π / 2.0) / coeff[1]
    return fn, fn¹

end

"""
    mathieu_Mce(kind, order, q, z)
    mathieu_Mce(kind, order, q, z::AbstractArray)

Computes the Even Modified Mathieu function `mathieu_Mce` and its derivative of `kind` (1 or 2) of order `m`, parameter `q` at `z`.
Computes the eigenvalue and coefficients only once when `z` is an array. 

"""
function mathieu_Mce(kind, order, q, z)
    (q <= 0) && @error "q must be greater than zero"
    ((kind < 1) || (kind > 2)) && @error "kind must be 1 or 2"
    # Compute the characteristic value and coefficients
    a = MathieuCharA(order, q)
    coeff = mathieu_a_coeff(order, q, a, 100)
    fn, fn¹ = Mce_kernel(kind, order, coeff, q, z; maxerr=1e-14)
    return fn, fn¹
end

function mathieu_Mce(kind, order, q, z::AbstractArray)
    (q <= 0) && @error "q must be > 0"
    ((kind < 1) || (kind > 2)) && @error "kind must be 1 or 2"
    # Compute the characteristic value and coefficients
    a = MathieuCharA(order, q)
    coeff = mathieu_a_coeff(order, q, a, 100)
    f = Mce_kernel.(kind, order, Ref(coeff), q, z; maxerr=1e-14)
    return f
end

"""
    MSe_kernel(kind, order, coeff, q, z; maxerr=1e-14)

Kernel for computing the odd Modified Mathieu function `mathieu_Mse` and its derivative of `kind` (1 or 2) of order `m`, parameter `q` at `z`.
"""
function Mse_kernel(kind, order, coeff, q, z; maxerr=1e-14)

    order == 0 && return (0.0, 0.0)
    amax = 0.0
    fn = 0.0
    fn¹ = zero(fn)
    u1 = sqrt(q) * exp(-z)
    u2 = sqrt(q) * exp(z)
    du1_dz = -u1
    du2_dz = u2

    if iseven(order)
        for k in 0:min(100, length(coeff) - 1)
            amax = max(amax, abs(coeff[k+1]))
            abs(coeff[k+1]) / amax < maxerr && break

            j1mc = besselj(k, u1)
            j1pc = besselj(k + 2, u1)
            j1mc_prime = besselj_prime(k, u1)
            j1pc_prime = besselj_prime(k + 2, u1)

            fl = kind == 1
            z2mc = fl ? besselj(k, u2) : bessely(k, u2)
            z2pc = fl ? besselj(k + 2, u2) : bessely(k + 2, u2)
            z2mc_prime = fl ? besselj_prime(k, u2) : bessely_prime(k, u2)
            z2pc_prime = fl ? besselj_prime(k + 2, u2) : bessely_prime(k + 2, u2)

            fc = (-1.0)^(0.5 * order + k + 1) * coeff[k+1]
            fn += fc * (j1mc * z2pc - j1pc * z2mc)
            fn¹ += fc * (j1mc_prime * z2pc * du1_dz + j1mc * z2pc_prime * du2_dz -
                         j1pc_prime * z2mc * du1_dz - j1pc * z2mc_prime * du2_dz)
        end
    else
        for k in 0:min(100, length(coeff) - 1)
            amax = max(amax, abs(coeff[k+1]))
            abs(coeff[k+1]) / amax < maxerr && break

            j1c = besselj(k, u1)
            j1pc = besselj(k + 1, u1)
            j1c_prime = besselj_prime(k, u1)
            j1pc_prime = besselj_prime(k + 1, u1)
            fl = kind == 1
            z2c = fl ? besselj(k, u2) : bessely(k, u2)
            z2pc = fl ? besselj(k + 1, u2) : bessely(k + 1, u2)
            z2c_prime = fl ? besselj_prime(k, u2) : bessely_prime(k, u2)
            z2pc_prime = fl ? besselj_prime(k + 1, u2) : bessely_prime(k + 1, u2)
            fc = (-1.0)^(0.5 * (order - 1) + k) * coeff[k+1]
            fn += fc * (j1c * z2pc - j1pc * z2c)
            fn¹ += fc * (j1c_prime * z2pc * du1_dz + j1c * z2pc_prime * du2_dz
                         -
                         j1pc_prime * z2c * du1_dz - j1pc * z2c_prime * du2_dz)
        end
    end
    fn *= sqrt(π / 2.0) / coeff[1]
    fn¹ *= sqrt(π / 2.0) / coeff[1]
    return fn, fn¹
end

"""
    mathieu_Mse(kind, order, q, z)
    mathieu_Mse(kind, order, q, z::AbstractArray)

Computes the Even Modified Mathieu function `mathieu_Mse` and its derivative of `kind` (1 or 2) of order `m`, parameter `q` at `z`.
Computes the eigenvalue and coefficients only once when `z` is an array. 

"""
function mathieu_Mse(kind, order, q, z)
    (q <= 0) && @error "q must be greather than zero"
    ((kind < 1) || (kind > 2)) && @error "kind mus tbe 1 or 2"
    # Compute the characteristic value
    a = MathieuCharB(order, q)
    coeff = mathieu_b_coeff(order, q, a, 100)
    fn, fn¹ = Mse_kernel(kind, order, coeff, q, z; maxerr=1e-14)
    return fn, fn¹
end

function mathieu_Mse(kind, order, q, z::AbstractArray)
    (q <= 0) && @error "q must be greather than zero"
    ((kind < 1) || (kind > 2)) && @error "kind mus tbe 1 or 2"
    # Compute the characteristic value
    a = MathieuCharB(order, q)
    coeff = mathieu_b_coeff(order, q, a, 100)
    f = Mse_kernel.(kind, order, Ref(coeff), q, z; maxerr=1e-14)
    return f
end


fc_even(n) = 4 * n^2
fc_odd(n) = (2 * n + 1)^2
fs_even(n) = 4 * (n + 1)^2
fs_odd(n) = (2 * n + 1)^2

function backward_recurse!(func::F, factors, a_param, q_param, xi, g, even_odd, parity, ni, N) where F
    factors[ni+1] = xi
    for n in 0:ni-1
        nn = N - n - 1
        factors[ni-n] = -1.0 / (((func(nn) - a_param)) / q_param + factors[ni-n+1])
    end
    if parity == false && even_odd == false && (ni == N - 1) # type: Mc/Ms, even_odd = even/odd m, ...
        factors[1] *= 2.0
    end
    return factors[1] - g
end

function backward_call!(fs_odd, factors, a_param, q_param, xi, g, even_odd, parity, ni, N)
    if parity == 0 # Mc / Ms
        if even_odd == 0 
            return backward_recurse!(fc_even, factors, a_param, q_param, xi, g, even_odd, parity, ni, N)
        else
            return backward_recurse!(fc_odd, factors, a_param, q_param, xi, g, even_odd, parity, ni, N)
        end
    else
        if even_odd == 0
            return backward_recurse!(fs_even, factors, a_param, q_param, xi, g, even_odd, parity, ni, N)
        else
            return backward_recurse!(fs_odd, factors, a_param, q_param, xi, g, even_odd, parity, ni, N)
        end
    end
end

function stabilize_backward_recurrence!(factors, a_param, q_param, even_odd, parity, ratio, ni, N; max_its = 20, eps = 1.0e-14)

    if parity == 0
        x1 = even_odd == 0 ? -q_param / fc_even(N) : -q_param / fc_odd(N)
    else
        x1 = even_odd == 0 ? -q_param / fs_even(N) : -q_param / fs_odd(N)
    end

    g1 = backward_call!(fs_odd, factors, a_param, q_param, x1, ratio, even_odd, parity, ni, N)

    x2 = g1
    g2 = backward_call!(fs_odd, factors, a_param, q_param, x1, ratio, even_odd, parity, ni, N)

    its = 0
    while its ≤ max_its

        error1 = g1 - x1
        error2 = g2 - x2
        Δerror = error1 - error2
       abs(Δerror) < eps && break

        xh = (error1 * x2 - error2 * x1) / Δerror

        x1, g1, x2, g2 = x2, g2, xh, ratio
        g2 = backward_call!(fs_odd, factors, a_param, q_param, x2, g2, even_odd, parity, ni, N)
        its += 1
    end

    return nothing
end


"""
    forward_recurrence!(coeff, order, a_param, q_param, even_odd, parity)

Compute the initial forward coefficients.

# Reference:

1. Algorithms for the Computation of all Mathieu Functions of Integer Orders, Fayez A. Alhargan.
2. Implementation of angular and radial Mathieu function based on GSL.

"""
function forward_recurrence!(coeff, order, a_param, q_param, even_odd, parity)
    if parity == 0
        return forward_recurrence_c!(coeff, order, a_param, q_param, even_odd)
    else
        return forward_recurrence_s!(coeff, order, a_param, q_param, even_odd)
    end
end


"""
    forward_recurrence_c!(coeff, order, a_param, q_param, even_odd)

Compute the initial forward coefficients for the even case (Mce)
"""
function forward_recurrence_c!(coeff, order, a_param, q_param, even_odd)

    sum_val = 0.0
    if even_odd == 0
        coeff[2] = a_param / q_param
        coeff[3] = (a_param - 4) / q_param * coeff[2] - 2.0
        sum_val = coeff[1] + coeff[2] + coeff[3]
        ii = 3
        while ii < order ÷ 2 + 1
            coeff[ii+1] = (a_param - 4 * (ii - 1)^2) / q_param * coeff[ii] - coeff[ii-1]
            sum_val += coeff[ii+1]
            ii += 1
        end
    else
        coeff[1+1] = (a_param - 1) / q_param - 1.0
        sum_val = coeff[0+1] + coeff[1+1]

        ii = 2
        while ii < order ÷ 2 + 1
            coeff[ii+1] = (a_param - (2 * ii - 1)^2) / q_param * coeff[ii] - coeff[ii-1]
            sum_val += coeff[ii+1]
            ii += 1
        end
    end

    nn = ii - 1
    ratio = coeff[nn+1] / coeff[nn]

    return ratio, nn, sum_val
end

"""
    forward_recurrence_s!(coeff, order, a_param, q_param, even_odd)

Compute the initial forward coefficients for the odd case (Mse)
"""
function forward_recurrence_s!(coeff, order, a_param, q_param, even_odd)

    sum_val = 0.0
    if even_odd == 0
        coeff[2] = (a_param - 4) / q_param
        sum_val += 2 * coeff[1] + 4 * coeff[2]
        ii = 2
        while ii < order ÷ 2
            coeff[ii+1] = (a_param - 4 * ii^2) / q_param * coeff[ii] - coeff[ii-1]
            sum_val += 2 * (ii + 1) * coeff[ii+1]
            ii += 1
        end
    else
        coeff[2] = (a_param - 1) / q_param + 1.0
        sum_val = coeff[1] + 3 * coeff[2]
        ii = 2
        while ii < order ÷ 2 + 1
            coeff[ii+1] = (a_param - (2 * ii - 1)^2) / q_param * coeff[ii] - coeff[ii-1]
            sum_val += (2 * (ii + 1) - 1) * coeff[ii+1]
            ii += 1
        end
    end

    nn = ii - 1
    ratio = coeff[nn+1] / coeff[nn]

    return ratio, nn, sum_val
end

"""
    normalize_coeffs!(coeff, factors, nn, sum_val, parity, N)

Extend coefficients with stabilized values and normalize.
"""
function normalize_coeffs!(coeff, factors, nn, sum_val, parity, N)

    if parity == 0
        sum_val += coeff[nn+1]
        for ii in (nn + 1):N-1
            coeff[ii+1] = factors[ii - nn] * coeff[ii]
            sum_val += coeff[ii+1]
            if abs(coeff[ii]) < 1e-20
                coeff[ii:end] .= 0.0
                break
            end
        end
    else
        sum_val += 2 * (nn + 1) * coeff[nn+1]
        for ii in (nn + 1):N-1
            coeff[ii+1] = factors[ii - nn] * coeff[ii]
            sum_val += 2 * (ii + 1) * coeff[ii+1]
            if abs(coeff[ii]) < 1e-20
                coeff[ii:end] .= 0.0
                break
            end
        end
    end
    coeff ./= sum_val
    return coeff
end

"""
    trivial_cases!(coeff, q_param, order, parity)

Handle trivial case for the computation of `a` and `b` coefficients.
"""
function trivial_cases!(coeff, q_param, order, parity)
    if parity == 0
        if q_param == 0.0
            coeff[order÷2+1] = 1.0
            return true
        end
    else
        if q_param == 0.0
            coeff[(order-1)÷2+1] = 1.0  #  1-based
            return true
        end
    end
    return false
end

"""
    low_order(a_param, q_param, even_odd, parity)

Handle low order cases for the computation of `a` and `b` coefficients.
"""
function low_order(a_param, q_param, even_odd, parity)
    nn = 0
    sum_val = 0.0
    if even_odd == 0
        ratio = parity == 0 ?  a_param / q_param : (a_param - 4) / q_param
    else
        ratio = (a_param - 1 - q_param) / q_param
    end
    return ratio, nn, sum_val
end

"""
    mathieu_a_coeff(order, q_param, a_param, N)

Computes Fourier coefficients for even Mathieu and Modified Mathieu functions.
"""
function mathieu_a_coeff(order, q_param, a_param, N)

    coeff = zeros(typeof(a_param), N,)
    coeff[1] = 1.0
    even_odd = isodd(order) ? 1 : 0
    #order > N && @error "Coefficient array is not large enough to hold all necessary coefficients"
    # Handle trivial cases
    trivial_cases!(coeff, q_param, order, 0) && return coeff

    # Forward recurrence
    if order < 5
        ratio, nn, sum_val = low_order(a_param, q_param, even_odd, false)
    else
        ratio, nn, sum_val = forward_recurrence!(coeff, order, a_param, q_param, even_odd, false)
    end

    ni = N - nn - 1
    # Backwards recurrence
    factors = zeros(typeof(a_param), N,)
    stabilize_backward_recurrence!(factors, a_param, q_param, even_odd, false, ratio, ni, N)
    
    normalize_coeffs!(coeff, factors, nn, sum_val, false, N)

    return coeff
end


"""
    mathieu_b_coeff(order, q_param, a_param, N)

Computes Fourier coefficients for Odd Mathieu and Modified Mathieu functions.
"""
function mathieu_b_coeff(order, q_param, a_param, N)

    coeff = zeros(typeof(a_param), N,)
    coeff[1] = 1.0
    even_odd = isodd(order) ? 1 : 0
    #order > N && @error "Coefficient array is not large enough to hold all necessary coefficients"
    # Handle trivial cases
    trivial_cases!(coeff, q_param, order, true) && return coeff

    # Forward recurrence
    if order < 5
        ratio, nn, sum_val = low_order(a_param, q_param, even_odd, true)
    else
        ratio, nn, sum_val = forward_recurrence!(coeff, order, a_param, q_param, even_odd, true)
    end
    ni = N - nn - 1
    # Backwards recurrence
    factors = zeros(typeof(a_param), N,)
    stabilize_backward_recurrence!(factors, a_param, q_param, even_odd, true, ratio, ni, N)
    # Compute the rest of the coefficients
    normalize_coeffs!(coeff, factors, nn, sum_val, true, N)

    return coeff
end

"""
    find_zero_ewg(func::F, m, q; inc_iters = 500, max_iters=70, tol=1e-10)

    Finds a numerical zero of `func` for a given initial estimate `q` and azimuthal index `m`.

This method is designed for eigenvalue problems in **elliptical waveguides (EWG)**, where directly using standard root-finding libraries (e.g. `Roots.jl`) can be unreliable without dense sampling or problem-specific knowledge.  
It combines a **sign-change search** and a **bisection refinement**, providing robust convergence for oscillatory or slowly decaying modal functions.

"""
function find_zero_ewg(func::F, m, q; inc_iters = 500, max_iters=70, tol=1e-10) where F
    factor = 1.1
    increment = 0.25
    x1 = q
    x2 = x1 * factor
    if isapprox(func(x1)[1], 0.0, atol = 1e-7)
        x1 = x1 * factor
        factor += increment
        x2 = x1 * factor
        f_a = func(x1)[1]
        f_b = func(x2)[1]
        sold = sign(f_a)
        snew = sign(f_b)
    else
        f_a = func(x1)[1]
        f_b = func(x2)[1]
        sold = sign(f_a)
        snew = sign(f_b)
    end

    # sign change
    its = 0
    while sold == snew && inc_iters > its
        sold = snew
        f_a = f_b
        x1 = x2
        factor += increment
        x2 = q * factor
        f_b = func(x2)[1]
        snew = sign(f_b)
        its += 1
    end

    # Bisection Method
    its = 0
    xc = (x1 + x2)/2
    while its <= max_iters
        xc = (x1 + x2)/2
        f_c = func(xc)[1]
        abs(f_c) <= tol && break
        if sign(f_a) == sign(f_c)
            x1 = xc
            f_a = f_c
        else
            x2 = xc
            f_b = f_c
        end
        its += 1
    end
    # Newton Method [No hace falta de momento]
    #its = 0
    #while its <= max_iters
    #    f, f¹ = func(xc)
    #    abs(f) <= tol && break
    #    xc -= f - f / f¹
    #    its += 1
    #end
    return xc
end

"""
    find_zero_ewg(func::F, m, T::Symbol; inc_iters = 500, max_iters = 70, tol = 1e-10)

Wrapper for [`find_zero_ewg`](@ref) that automatically selects an initial guess `q`
based on the mode type `T` (`:TE` or `:TM`) and the azimuthal index `m`.
"""
function find_zero_ewg(func::F, m, T::Symbol; inc_iters = 500, max_iters=70, tol=1e-10) where F
    if T == :TE
        q = max(0.1, m/2 * 0.7)
    else
        q = max(0.1, 1.8*m)
    end
    return find_zero_ewg(func, m, q; inc_iters = inc_iters, max_iters=max_iters, tol=tol)
end


"""
    find_n_zero_ewg(func::F, m, n, T; inc_iters = 500, max_iters = 70, tol = 1e-10)

Finds the `n`-th zero of `func` for a given azimuthal index `m` and mode type `T`.

Internally calls [`find_zero_ewg`](@ref) iteratively, using each previously found
zero as the starting point for the next one.
"""
function find_n_zero_ewg(func::F, m, n, T; inc_iters = 500, max_iters=70, tol=1e-10) where F
    n == 0 && return 0.0
    v = find_zero_ewg(func, m, T; inc_iters = inc_iters, max_iters=max_iters, tol=tol)
    for _ in 2:n
        v = find_zero_ewg(func, m, v; inc_iters = inc_iters, max_iters=max_iters, tol=tol)
    end
    return v
end
