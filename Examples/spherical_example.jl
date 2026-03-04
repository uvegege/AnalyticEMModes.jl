using FastGaussQuadrature, StaticArrays, LinearAlgebra

function sample_sphere(r0, Nθ, Nφ)
    # Gauss-Legendre nodes for θ in [0,π]
    x, w = gausslegendre(Nθ)             # x in [-1,1]
    thetas = acos.(x)                    # θ nodes
    wθ = w                               # weights correspond to dx = d(cosθ)
    phis = range(0, 2π, length=Nφ+1)[1:end-1]
    # build list of (θ, φ, weight)
    pts = SVector{3, Float64}[]
    for (i,θ) in enumerate(thetas)
        for φ in phis
            weight = wθ[i] * (2π / Nφ)  
            push!(pts, SVector(θ, φ, weight))
        end
    end
    return pts 
end

function sph2cart_field(E_r, E_θ, E_ϕ, θ, ϕ)
    st, ct = sincos(θ)
    sp, cp = sincos(ϕ)
    E_x = E_r * st*cp + E_θ * ct*cp - E_ϕ * sp
    E_y = E_r * st*sp + E_θ * ct*sp + E_ϕ * cp
    E_z = E_r * ct     - E_θ * st
    return SVector(E_x, E_y, E_z)
end

function cart2sph_field(E_x, E_y, E_z, θ, ϕ)
    st, ct = sincos(θ)
    sp, cp = sincos(ϕ)
    E_r = E_x * st * cp + E_y * st * sp + E_z * ct
    E_θ = E_x * ct * cp + E_y * ct * sp - E_z * st
    E_ϕ = -E_x * sp + E_y * cp
    return SVector(E_r, E_θ, E_ϕ)
end


# ------------ assemble modal matrix ----------
function assemble_modal_matrix(Lmax, k, r0, Nϕ, Nθ)

    points_sph = sample_sphere(r0, Nθ, Nϕ)
    points_cart = map(x->SVector(AnalyticEMModes.sph2cart(r0, x[1], x[2])), points_sph)
    weights = map(x->x[3], points_sph)
    modes = eachcol(mn_sph_vectors_lmax(points_cart, Lmax, 10e9, 1.0, 1.0, true)) # Each col is a mode (l,m) mode

    nmodes = 2*length(modes)
    P = length(points_cart)

    A = zeros(ComplexF64, 3*P, nmodes)
    b = zeros(ComplexF64, 3*P)

    # loop points
    row = 1
    for (pidx, (x, y, z)) in enumerate(points_cart)
        
        r, θ, ϕ = AnalyticEMModes.cart2sph(x, y, z)
        # evaluate test field on sphere (plane wave in +z, pol x)
        #Et = SVector(1.0 + 0im, 0.0 + 0im, 0.0 + 0im) * exp(im * k * z)  # plane wave, x-pol
        Et = SVector(sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))* exp(im * k * r)
        # store b with sqrt(weight) applied
        s = sqrt(weights[pidx])  # quadrature weight factor (no r^2 because we used d(cosθ)*dφ)
        b[row]   = s * Et[1]
        b[row+1] = s * Et[2]
        b[row+2] = s * Et[3]

        for (j, md) in enumerate(modes)
            mr, mθ, mϕ, nr, nθ, nϕ = md[pidx]
            mx, my, mz, nx, ny, nz = spherical_to_cartesian_fields(mr, mθ, mϕ, nr, nθ, nϕ, θ, ϕ)
            # put weighted mode components into A
            A[row, j]   = s * mx
            A[row+1, j] = s * my
            A[row+2, j] = s * mz
            A[row, j+length(modes)]   = s * nx
            A[row+1, j+length(modes)] = s * ny
            A[row+2, j+length(modes)] = s * nz
        end

        row += 3
    end

    return A, b, modes
end


function assemble_modal_matrix(Lmax, k, r0, Nϕ, Nθ)

    points_sph = sample_sphere(r0, Nθ, Nϕ)
    points_cart = map(x->SVector(AnalyticEMModes.sph2cart(r0, x[1], x[2])), points_sph)
    weights = map(x->x[3], points_sph)
    modes = eachcol(mn_sph_vectors_lmax(points_cart, Lmax, 10e9, 1.0, 1.0, true)) # Each col is a mode (l,m) mode

    nmodes = length(modes)
    P = length(points_cart)

    G = zeros(ComplexF64, 3*P, 2*nmodes)
    f = zeros(ComplexF64, 3*P)
    
    for i in 1:P
        s = sqrt(weights[i])
        r, θ, ϕ = points_sph[i]
        er, et, ep = cart2sph_field(cis(-k * r*cos(θ)), 0.0, 0.0, θ, ϕ)
        Et = SVector(0.0, cos(3*θ)^2, sin(3*θ)^2) * cis(-k * r*cos(θ))
        #f[3*i-2] = s * er
        #f[3*i-1] = s * et
        #f[3*i]   = s * ep
        #f[3*i-2] = s * (modes[4][i][1] + 1/3 * modes[15][i][1] - 1/8 * modes[27][i][4] )
        #f[3*i-1] = s * (modes[4][i][2] + 1/3 * modes[15][i][2] - 1/8 * modes[27][i][5] )
        #f[3*i]   = s * (modes[4][i][3] + 1/3 * modes[15][i][3] - 1/8 * modes[27][i][6] )
        f[3*i-2] = s * (modes[4][i][1]^2 + 1/3 * modes[15][i][1]^3)
        f[3*i-1] = s * (modes[4][i][2]^2 + 1/3 * modes[15][i][2]^3)
        f[3*i]   = s * (modes[4][i][3]^2 + 1/3 * modes[15][i][3]^3)
    end

    for (j, md) in enumerate(modes)
        for i in 1:P
            s = sqrt(weights[i])
            G[3*i - 2, j] = s * md[i][1] # Er/Mr/Nr
            G[3*i - 1, j] = s * md[i][2] # Eθ/Mθ/Nθ
            G[3*i, j]     = s * md[i][3] # Eϕ/Mϕ/Nϕ
            G[3*i - 2, j + length(modes)] = s * md[i][4] # Er/Mr/Nr
            G[3*i - 1, j + length(modes)] = s * md[i][5] # Eθ/Mθ/Nθ
            G[3*i, j + length(modes)]     = s * md[i][6] # Eϕ/Mϕ/Nϕ
        end
    end
   
    return G, f
end


function modal_decomposition(Lmax, f, r0; Nθ=81, Nφ=81)
    A, b = assemble_modal_matrix(Lmax, wavenumber(f, 1, 1), r0, Nθ, Nφ)
    c = A \ b  
    res = norm(A*c - b) / norm(b)   # relative residual
    return c, res
end

modal_decomposition(21, 10e9, 1e-3; Nθ=81, Nφ=81)


function inner_product(mode1, mode2, points_sph, weights)
    S = 0.0 + 0im
    for i in eachindex(points_sph)
        _, θ, ϕ = points_sph[i]
        mr1, mth1, mph1 = mode1[i]
        mr2, mth2, mph2 = mode2[i]
        #(mx1,my1,mz1) = sph2cart_field(mr1, mth1, mph1, θ, ϕ)
        #(mx2,my2,mz2) = sph2cart_field(mr2, mth2, mph2, θ, ϕ)
        (mx1,my1,mz1) = mr1, mth1, mph1
        (mx2,my2,mz2) = mr2, mth2, mph2
        w = weights[i]
        S += w * (mx1*conj(mx2) + my1*conj(my2) + mz1*conj(mz2))
    end
    return S
end

function inner_product_test(mode1, mode2, points_sph, weights)
    S = 0.0 + 0im
    for i in eachindex(points_sph)
        _, θ, ϕ = points_sph[i]
        mr1, mth1, mph1 = mode1[i]
        mr2, mth2, mph2 = mode2[i]
        #(mx1,my1,mz1) = sph2cart_field(mr1, mth1, mph1, θ, ϕ)
        #(mx2,my2,mz2) = sph2cart_field(mr2, mth2, mph2, θ, ϕ)
        (mx1,my1,mz1) = mr1, mth1, mph1
        (mx2,my2,mz2) = mr2, mth2, mph2
        w = weights[i]
        S += w * (mx1*conj(mx2) + my1*conj(my2) + mz1*conj(mz2))
    end
    return S
end


begin
    Nθ=51
    Nφ=51
    r0 = 12e-3
    freq = 10.0e9
    Lmax = 7
    k = wavenumber(freq, 1, 1)
    points_sph = sample_sphere(r0, Nθ, Nφ)
    points_cart = map(x->SVector(AnalyticEMModes.sph2cart(r0, x[1], x[2])), points_sph)
    weights = map(x->x[3], points_sph)
    modes = eachcol(mn_sph_vectors_lmax(points_cart, Lmax, freq, 1.0, 1.0, true));
    Mmodes = map(x->map(y->SVector(y[1], y[2], y[3]), x), modes);
    Nmodes = map(x->map(y->SVector(y[4], y[5], y[6]), x), modes);

    A = zeros(ComplexF64, length(modes), length(modes))
    weights = map(x -> x[3], points_sph)
    for i in eachindex(modes)
        A[i,i] = inner_product(Mmodes[i], Mmodes[i], points_sph, weights)
    end

    maximum(abs.(A - Diagonal(diag(A))))
    abs.(Diagonal(A).diag)
end

spy(abs.(A))

begin
    A = zeros(ComplexF64, length(modes), length(modes))
    weights = map(x -> x[3], points_sph)
    for i in eachindex(modes)
        for j in eachindex(modes)
            A[i,j] = inner_product(Mmodes[i], Mmodes[j], points_sph, weights)
        end
    end
    maximum(abs.(A - Diagonal(diag(A))))
end

spy(abs.(A))

begin
    A = zeros(ComplexF64, length(modes), length(modes))
    weights = map(x -> x[3], points_sph)
    for i in eachindex(modes)
        for j in eachindex(modes)
            A[i,j] = inner_product(Mmodes[i], Nmodes[j], points_sph, weights)
        end
    end
    maximum(abs.(A))
end


spy(abs.(A))

begin
    A = zeros(ComplexF64, length(modes), length(modes))
    weights = map(x -> x[3], points_sph)
    for i in eachindex(modes)
        for j in eachindex(modes)
            A[i,j] = inner_product(Nmodes[i], Nmodes[j], points_sph, weights)
        end
    end
    maximum(abs.(A))
end
spy(abs.(A))


# ------------ assemble modal matrix ----------


function test1(points)
    return [SVector(1.0 + 0im, 0.0 + 0im, 0.0 + 0im) * exp(im * k * z) for (x,y,z) in points]
end

function test1(points)
    return [SVector(sin(y) * cos(z),sin(y) * sin(z), cos(y))  for (x,y,z) in points]
end

function assemble_modal_matrix_test(Lmax, r0, Nϕ, Nθ)

    points_sph = sample_sphere(r0, Nθ, Nϕ)
    points_cart = map(x->SVector(AnalyticEMModes.sph2cart(r0, x[1], x[2])), points_sph)
    weights = map(x->x[3], points_sph)
    modes = eachcol(mn_sph_vectors_lmax(points_cart, Lmax, 10e9, 1.0, 1.0, true)) # Each col is a mode (l,m) mode
    Mmodes = map(x->map(y->SVector(y[1], y[2], y[3]), x), modes)
    Nmodes = map(x->map(y->SVector(y[4], y[5], y[6]), x), modes)
    sol = test1(points_sph)
    #sol = map(eachindex(points_sph)) do i
    #    return SVector((modes[4][i][1] + 1/3 * modes[15][i][1] - 1/8 * modes[27][i][4] ),
    #    (modes[4][i][2] + 1/3 * modes[15][i][2] - 1/8 * modes[27][i][5] ),
    #    (modes[4][i][3] + 1/3 * modes[15][i][3] - 1/8 * modes[27][i][6] ))
    #end

    #sol = map(eachindex(points_sph)) do i
    #    return SVector((modes[4][i][1] + 1/3 * modes[15][i][1]),
    #    (modes[4][i][2] + 1/3 * modes[15][i][2]),
    #    (modes[4][i][3] + 1/3 * modes[15][i][3]))
    #end

    Alm = zeros(ComplexF64, length(Mmodes), )
    Blm = zeros(ComplexF64, length(Nmodes), )
    for i in eachindex(modes)
        Alm[i] = inner_product_test(Mmodes[i], sol, points_sph, weights)
    end
    for i in eachindex(modes)
        Blm[i] = inner_product_test(Nmodes[i], sol, points_sph, weights)
    end

    return Alm, Blm, modes, sol
end

begin
    Alm, Blm, modes, sol = assemble_modal_matrix_test(21, 1e-3, 81, 81)
    sol_modes = zeros(MVector{3, ComplexF64}, size(sol))
    Mmodes = map(x->map(y->SVector(y[1], y[2], y[3]), x), modes)
    Nmodes = map(x->map(y->SVector(y[4], y[5], y[6]), x), modes)
    sol_modes = sum(Alm .* Mmodes) #+ sum(Blm .* Nmodes)
    #@show norm(sol_modes + sol)
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, abs.(getindex.(sol, 1)))
    lines!(ax, abs.(getindex.(sol_modes, 1)))
    ax = Axis(fig[1,2])
    lines!(ax, abs.(getindex.(sol, 2)))
    lines!(ax, abs.(getindex.(sol_modes, 2)))
    ax = Axis(fig[1,3])
    lines!(ax, abs.(getindex.(sol, 3)))
    lines!(ax, abs.(getindex.(sol_modes, 3)))
    fig
end
