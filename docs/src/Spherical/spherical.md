```@meta
    ShareDefaultModule = true
```

```@setup
include("../Assets/docs_include.jl")
```


## TE Modes

For the computation of spherical modes, for reasons of efficiency and implementation, it is recommended to use a vector `Vector{SVector{3,T}}` and calculate `(lmax+1)^2` modes directly. 

With `lmax = 6` you have 49 **TE** modes.

```@example

lmax = 6
R = 1
name = "../Assets/mesh/sphregion.msh"
#sphregion_mesh(R, dl = 0.75e-1, name = name)

coord, conn = mesh_data(name)
coords = coord[:, 1:maximum(conn)]
xcoords = coords[1, :]
ycoords = coords[2, :]
zcoords = coords[3, :]
#r_vec = [SVector(v...) for v in eachcol(coords)]
r_vec = to_svector(xcoords, ycoords, zcoords)


myfields = te_sph_fields_lmax(r_vec, lmax, 10e9, 1.0, 1.0, 4)

lm_pairs = [(4, 1), (6, -2)]

plot_ids = Iterators.product(1:2, 1:1) |> collect |> vec

fig = Figure(size = (900, 800))
for (idj, (l, m)) in enumerate(lm_pairs)
    ii, jj = plot_ids[idj]
    #l, m = 3, 2
    idx = l^2 + l + m + 1

    mode = myfields[:,idx]
    E = map(mode, r_vec) do mi, ri
        E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ = mi
        r, θ, ϕ = AnalyticEMModes.cart2sph(ri...)
        spherical_to_cartesian_fields(0.0, E_θ, E_ϕ, 0.0, H_θ, H_ϕ, θ, ϕ)
    end

    Hx = map(x->real(x[4]), E)
    Hy = map(x->real(x[5]), E)
    Hz = map(x->real(x[6]), E)
    Hf = (map((x,y,z) -> hypot(x,y,z), Hx, Hy, Hz))
    Hx ./= Hf
    Hy ./= Hf
    Hz ./= Hf
    factor = 0.1

ga = fig[ii, jj] = GridLayout()
    Label(ga[0, 1:2], L"TE_{%$l,%$m}", fontsize = 30, tellwidth = false, padding = (0, 0, -170, 0))
    Label(ga[0, 1], L"H_t", fontsize = 30, tellwidth = false, padding = (0, 0, -170, 0))
    ax1 = LScene(ga[1,1], show_axis = false)
    mesh!(ax1, coords, conn, color = :gray)
    arrows3d!(ax1, xcoords, ycoords, zcoords, Hx*factor, Hy*factor, Hz*factor, color = Hf, colormap = :jet)
    Label(ga[0,2], L"H_\rho", fontsize = 30, tellwidth = false, padding = (0, 0, -170, 0))
    ax2 = LScene(ga[1,2], show_axis = false)
    Hρ = map(x->real(x[4]), mode)
    mesh!(ax2, coords, conn, color = real.(Hρ), colormap = :jet)
    rotate_cam!(ax2.scene, Vec3f(0, 0.0, 0.0))
    rotate_cam!(ax1.scene, Vec3f(0, 0.0, 0.0))
end
fig
```

## M/N Basis Visualization

You can also visualize the raw spherical vector-wave basis directly.

```@example
lmax = 4
f = 10e9
μr, εr = 1.0, 1.0
radial = 4
name = "../Assets/mesh/sphregion.msh"

coord, conn = mesh_data(name)
coords = coord[:, 1:maximum(conn)]
xcoords = coords[1, :]
ycoords = coords[2, :]
zcoords = coords[3, :]
r_vec = to_svector(xcoords, ycoords, zcoords)

mn = mn_sph_vectors_lmax(r_vec, lmax, f, μr, εr, radial)

l, m = 3, 1
idx = l^2 + l + m + 1
mode = mn[:, idx]

Mcart = map(mode, r_vec) do mi, ri
    Mr, Mθ, Mϕ, Nr, Nθ, Nϕ = mi
    _, θ, ϕ = AnalyticEMModes.cart2sph(ri...)
    Ex, Ey, Ez, _, _, _ = spherical_to_cartesian_fields(Mr, Mθ, Mϕ, 0.0, 0.0, 0.0, θ, ϕ)
    (Ex, Ey, Ez, Nr)
end

Mx = map(v -> real(v[1]), Mcart)
My = map(v -> real(v[2]), Mcart)
Mz = map(v -> real(v[3]), Mcart)
Nr = map(v -> real(v[4]), Mcart)
Mf = map((x, y, z) -> hypot(x, y, z), Mx, My, Mz)

Mx ./= Mf
My ./= Mf
Mz ./= Mf

fig = Figure(size = (900, 420))
ga = fig[1, 1] = GridLayout()

Label(ga[0, 1:2], L"M_{%$l,%$m}\ \mathrm{and}\ N_r", fontsize = 28, tellwidth = false)
ax1 = LScene(ga[1, 1], show_axis = false)
mesh!(ax1, coords, conn, color = :gray)
arrows3d!(ax1, xcoords, ycoords, zcoords, 0.08 .* Mx, 0.08 .* My, 0.08 .* Mz, color = Mf, colormap = :turbo)

ax2 = LScene(ga[1, 2], show_axis = false)
mesh!(ax2, coords, conn, color = Nr, colormap = :balance)

rotate_cam!(ax1.scene, Vec3f(0.0, 0.0, 0.0))
rotate_cam!(ax2.scene, Vec3f(0.0, 0.0, 0.0))
fig
```
