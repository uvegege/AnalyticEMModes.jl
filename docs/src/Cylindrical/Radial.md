```@meta
    ShareDefaultModule = true
```

```@setup
include("../Assets/docs_include.jl")
```

# Radial Waveguide

Radial waveguides are structures that support radial wave propagation between parallel conducting plates. Unlike conventional waveguides where propagation is along the z-axis, radial guides exhibit wave propagation in the radial direction (ρ-direction). The field solutions involve Hankel functions of the radial coordinate, combined with sinusoidal variations in the azimuthal (φ) and vertical (z) directions.

Common examples include parallel-plate disk structures and sector-shaped geometries. For detailed mathematical derivations, see the Theory section.

## TE Modes

The following visualization shows the Hρ field for the first 12 TE modes, rendered on a 3D cylindrical sector:

```@example

H = 1
R = 1

name = "../Assets/mesh/radial_wg1.msh"
#radialwg_mesh(H, R; dl = 0.4e-1, name = name)

coord, conn, _ = mesh_data_radial(name)
coords = coord[:, 1:maximum(conn)]
xcoords = coords[1, :]
ycoords = coords[2, :]
zcoords = coords[3, :]
rcoords = map((x, y) -> hypot(x, y), xcoords, ycoords)
ϕcoords = map((x, y) -> atan(y, x), xcoords, ycoords)

modekind = [(0, 1), (1, 1), (2, 1), (3, 1), (0, 2), (1, 2),
            (2, 2), (3, 2), (0, 3), (1, 3), (2, 3), (3, 3)]

plot_ids = Iterators.product(1:4, 1:3) |> collect |> vec

fig = Figure(size = (1200, 900))
for (idplot, (m, n)) in enumerate(modekind)
    stitle = L"TE_{%$m%$n}"
    ii, jj = plot_ids[idplot]
    axi = Axis3(fig[jj,ii], title = stitle, titlesize = 20)
    hidedecorations!(axi)
    fields = te_radial_fields(rcoords, ϕcoords, zcoords, H, m, n, 0.0, 1.0, 100e9, 1, 1)
    fz = getindex.(fields, 4)
    mesh!(axi, coords, conn, color = abs.(fz), colormap = :jet, interpolate = false)
end
fig
```

The mode indices (m, n) determine the number of azimuthal and vertical variations, respectively. Notice how modes with m = 0 are axisymmetric, while modes with m > 0 exhibit azimuthal nodal planes.

## TM Modes

The following shows the Eρ field patterns for the first 12 TM modes:

```@example

fig = Figure(size = (1200, 900))

modekind = [(0, 1), (1, 1), (2, 1), (3, 1), 
            (0, 2), (1, 2), (2, 2), (3, 2), 
            (0, 3), (1, 3), (2, 3), (3, 3)]

plot_ids = Iterators.product(1:4, 1:3) |> collect |> vec
for (idplot, (m, n)) in enumerate(modekind)
    stitle = L"TM_{%$m%$n}"
    ii, jj = plot_ids[idplot]
    axi = Axis3(fig[jj,ii], title = stitle, titlesize = 20)
    hidedecorations!(axi)
    fields = tm_radial_fields(rcoords, ϕcoords, zcoords, H, m, n, 0.0, 1.0, 100e9, 1, 1)
    fz = getindex.(fields, 1)
    mesh!(axi, coords, conn, color = abs.(fz), colormap = :jet, interpolate = false)
end
fig
```

TM modes satisfy boundary conditions where Eρ = 0 at the top and bottom conducting plates (at z = 0 and z = H).

# Wedge Waveguide

Wedge waveguides are angular sectors bounded by two radial conducting walls at angles φ = 0 and φ = φ₀, along with top and bottom parallel plates. These structures support modes with quantized azimuthal variations determined by the wedge angle. The mode parameter m must satisfy m = nπ/φ₀ where n is an integer, ensuring proper boundary conditions at the angular walls.

The following visualization shows the Hz field (longitudinal magnetic component) for the first 12 TE modes in a wedge with angle φ₀ = π/3 (60 degrees):

```@example

H = 1
R = 1

name = "../Assets/mesh/wedge_wg1.msh"

ϕ0 = pi/3
#wedgewg_mesh(H, R, ϕ0; dl = 0.3e-1, name = name)
coord, conn, connf = mesh_data_radial(name)
c_ids = unique(conn)
coords = coord[:, c_ids]
D = Dict(c_ids .=> eachindex(c_ids))
new_conn = map(x->D[x], conn)
xcoords = coords[1, :]
ycoords = coords[2, :]
zcoords = coords[3, :]
rcoords = map((x, y) -> hypot(x, y), xcoords, ycoords)
ϕcoords = map((x, y) -> atan(y, x), xcoords, ycoords)


modekind = [(0, 1), (1, 1), (2, 1), (3, 1), 
            (0, 2), (1, 2), (2, 2), (3, 2), 
            (0, 3), (1, 3), (2, 3), (3, 3)]

plot_ids = Iterators.product(1:4, 1:3) |> collect |> vec

fig = Figure(size = (1200, 900))
for (idplot, (m, n)) in enumerate(modekind)
    stitle = L"TE_{%$m%$n}"
    ii, jj = plot_ids[idplot]
    axi = Axis3(fig[jj,ii], title = stitle, titlesize = 20, azimuth = pi/5)
    hidedecorations!(axi)
    fields = te_wedge_fields(rcoords, ϕcoords, zcoords, H, ϕ0, m, n, 1.0, 0.0, 100e9, 1, 1)
    fz = getindex.(fields, 6)
    m1 = mesh!(axi, coord, connf, color = :gray)
    m2 = mesh!(axi, coords, new_conn, color = abs.(fz), colormap = :jet)
end
fig
```

