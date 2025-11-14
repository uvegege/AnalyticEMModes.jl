function rwg_mesh(W, H; dl = 1e-1, name = "./rectangular_wg.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("rw")
    rw = gmsh.model.occ.add_rectangle(0, 0, 0, W, H)
    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end

    PEC = Gmsh.gmsh.model.getEntities(1)

    Gmsh.gmsh.model.addPhysicalGroup(1, getindex.(PEC, 2), 1, "PEC")
    
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(Gmsh.gmsh.model.getEntities(2), 2), 1, "Domain")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.write(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap only
    Gmsh.gmsh.finalize()
    return model
end


# Circular mesh


function cwg_mesh(R; dl = 1e-1, name = "./circular_wg.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("cw")
    cw = gmsh.model.occ.add_disk(0, 0, 0, R, R)
    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end

    PEC = Gmsh.gmsh.model.getEntities(1)

    Gmsh.gmsh.model.addPhysicalGroup(1, getindex.(PEC, 2), 1, "PEC")
    
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(Gmsh.gmsh.model.getEntities(2), 2), 1, "Domain")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.write(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap only
    Gmsh.gmsh.finalize()
    return model
end

# Coaxial mesh


function coaxwg_mesh(R1, R2; dl = 1e-1, name = "./coaxial_wg.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("coaxw")
    # = gmsh.model.occ.add_disk(0, 0, 0, R, R)
    l1 = gmsh.model.occ.add_circle(0,0,0,R1)
    cl1 = gmsh.model.occ.add_curve_loop([l1])
    l2 = gmsh.model.occ.add_circle(0,0,0,R2)
    cl2 = gmsh.model.occ.add_curve_loop([l2])
    coaxw = gmsh.model.occ.addPlaneSurface([cl1, cl2])

    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end

    PEC = Gmsh.gmsh.model.getEntities(1)

    Gmsh.gmsh.model.addPhysicalGroup(1, getindex.(PEC, 2), 1, "PEC")
    
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(Gmsh.gmsh.model.getEntities(2), 2), 1, "Domain")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.write(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap only
    Gmsh.gmsh.finalize()
    #return model
end


# Radial mesh


function radialwg_mesh(H, R; dl = 1e-1, name = "./radial_wg.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("radxw")
    radwg = gmsh.model.occ.add_cylinder(0, 0, 0, 0, 0, H, R)

    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end

    tags2D = Gmsh.gmsh.model.getEntities(2)
    Domain = filter(x->gmsh.model.getType(x...) == "Cylinder", tags2D)
    PEC_Surfaces = filter(x->gmsh.model.getType(x...) == "Plane", tags2D)
    
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(Domain, 2), 1, "Domain")
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(PEC_Surfaces, 2), 2, "PEC")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    Gmsh.gmsh.write(name)
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap Only
    Gmsh.gmsh.finalize()
    #return model
end

# Wedge mesh


function wedgewg_mesh(H, R, angle; dl = 1e-1, name = "./wedge_wg.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("wedge")

    p1 = gmsh.model.occ.add_point(R * cos(-angle/2), R * sin(-angle/2), H)
    p2 = gmsh.model.occ.add_point(R * cos(0), R * sin(0), H)
    p3 = gmsh.model.occ.add_point(R * cos(+angle/2), R * sin(+angle/2), H)

    p4 = gmsh.model.occ.add_point(R * cos(-angle/2), R * sin(-angle/2), 0)
    p5 = gmsh.model.occ.add_point(R * cos(0), R * sin(0), 0)
    p6 = gmsh.model.occ.add_point(R * cos(+angle/2), R * sin(+angle/2), 0)

    c1 = gmsh.model.occ.add_point(0.0, 0.0, H)
    c2 = gmsh.model.occ.add_point(0.0, 0.0, 0)

    arc1 = gmsh.model.occ.add_circle_arc(p1, p2, p3, -1, false)
    arc2 = gmsh.model.occ.add_circle_arc(p4, p5, p6, -1, false)

    vl1 = gmsh.model.occ.add_line(p1, p4)
    vl2 = gmsh.model.occ.add_line(p3, p6)
    vl3 = gmsh.model.occ.add_line(p1, c1)
    vl4 = gmsh.model.occ.add_line(p3, c1)
    vl5 = gmsh.model.occ.add_line(p4, c2)
    vl6 = gmsh.model.occ.add_line(p6, c2)
    vl7 = gmsh.model.occ.add_line(c1, c2)

    domain_face = gmsh.model.occ.addCurveLoop([arc1, vl1, -arc2, vl2])
    pec_face_top = gmsh.model.occ.addCurveLoop([arc1, vl4, -vl3])
    pec_face_bot = gmsh.model.occ.addCurveLoop([arc2, vl6, -vl5])
    pec_face_left = gmsh.model.occ.addCurveLoop([vl7, -vl5, -vl1, vl3])
    pec_face_right = gmsh.model.occ.addCurveLoop([vl7, -vl6, -vl2, vl4])

    top_face = gmsh.model.occ.addPlaneSurface([pec_face_top])
    bot_face = gmsh.model.occ.addPlaneSurface([pec_face_bot])
    left_face = gmsh.model.occ.addPlaneSurface([pec_face_left])
    right_face = gmsh.model.occ.addPlaneSurface([pec_face_right])
    radial_face = gmsh.model.occ.addSurfaceFilling(domain_face) 

    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end
 
    Gmsh.gmsh.model.addPhysicalGroup(2, [radial_face], 1, "Domain")
    Gmsh.gmsh.model.addPhysicalGroup(2, [top_face, bot_face, left_face, right_face], 2, "PEC")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    Gmsh.gmsh.write(name)
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap Only
    Gmsh.gmsh.finalize()
    #return model
end


# Elliptic mesh

function ellipticwg_mesh(R1, R2; dl = 1e-1, name = "./elliptic_wg.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("ew")
    ew = Gmsh.gmsh.model.occ.add_disk(0.0, 0.0, 0, R1, R2)
    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end

    PEC = Gmsh.gmsh.model.getEntities(1)

    Gmsh.gmsh.model.addPhysicalGroup(1, getindex.(PEC, 2), 1, "PEC")
    
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(Gmsh.gmsh.model.getEntities(2), 2), 1, "Domain")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.model.mesh.set_order(1)
    Gmsh.gmsh.write(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap only
    Gmsh.gmsh.finalize()
    #return model
end





function sphregion_mesh(R1; dl = 1e-1, name = "./sphregion.msh")

    Gmsh.gmsh.initialize()
    Gmsh.gmsh.option.setNumber("General.Terminal", 1)
    Gmsh.gmsh.clear()
    Gmsh.gmsh.model.add("ew")
    sph = Gmsh.gmsh.model.occ.add_sphere(0.0, 0.0, 0.0, R1)
    Gmsh.gmsh.model.occ.synchronize()

    for (dim, tag) in Gmsh.gmsh.model.getEntities(0)
        Gmsh.gmsh.model.mesh.setSize((dim, tag), dl)
    end

    
    Gmsh.gmsh.model.addPhysicalGroup(2, getindex.(Gmsh.gmsh.model.getEntities(2), 2), 1, "Domain")

    Gmsh.gmsh.model.occ.synchronize()
    Gmsh.gmsh.model.mesh.generate(2)
    Gmsh.gmsh.model.mesh.set_order(1)
    Gmsh.gmsh.write(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()
    #model = GmshDiscreteModel(Gmsh.gmsh) # Gridap only
    Gmsh.gmsh.finalize()
    #return model
end

sphregion_mesh(1)
coord, conn = mesh_data("sphregion.msh") 


lmax = 6
@btime te_sph_fields_lmax($r_vec, 6, 10e9, 1.0, 1.0);
@btime basis = SphericalHarmonics($lmax)
basis = SphericalHarmonics(lmax)
@btime Ylm, ∇Ylm = compute_with_gradients($basis, $r_vec)


@profview_allocs compute_with_gradients(basis, r_vec) [sample_rate = 0.1]


@profview_allocs sample_rate=1 compute_with_gradients(basis, r_vec)
@btime compute_with_gradients($basis, $r_vec)

@profview_allocs sample_rate=1 te_sph_fields_lmax(r_vec, 6, 10e9, 1.0, 1.0);


@profview te_sph_fields_lmax(r_vec, 6, 10e9, 1.0, 1.0)

@profview_allocs sample_rate=1 begin 
    allocfields = te_sph_fields_lmax(r_vec, 6, 10e9, 1.0, 1.0);
end


@code_warntype te_sph_fields_lmax(r_vec, 6, 10e9, 1.0, 1.0; incident = false)


prof = Profile.Allocs.fetch();


@profview_allocs begin
    for _ in 1:10
    te_sph_fields_lmax(r_vec, 6, 10e9, 1.0, 1.0)
    end
end

ijk = 2
l = 0
idx = 3
te_sph_fields(r_vec[ijk], Rs[ijk, l+1], Ylm[ijk, idx], ∇Ylm[ijk, idx], k, μᵣ, εᵣ)

coords = coord[:, 1:maximum(conn)]
xcoords = coords[1, :]
ycoords = coords[2, :]
zcoords = coords[3, :]
r_vec = [SVector(v...) for v in eachcol(coords)]
myfields = te_sph_fields_lmax(r_vec, 6, 10e9, 1.0, 1.0)
mode1 = myfields[:, 2]


sph_h1m_with_derivatives(0, 5*1.0) |> typeof

f = 10e9
μᵣ, εᵣ = 1.0, 1.0
lmax = 6
incident = false
A = Array{NTuple{6, ComplexF64}, 2}(undef, length(r_vec), (lmax+1)^2)
Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)


function sph_bessel_with_derivatives(f::F, l, x) where F
    h  = f(l, x)
    h′ = (l/x) * h - f(l+1, x)  
    #h″ = ((l*(l+1))/x^2 - 1) * h - (2/x) * h′
    # return (R``+k^2)Fr
    h″ = l*(l+1)/x^2 * h
    return (h, h′, h″)
end


jm(m, x) = sphericalbesselj(m, x)
ym(m, x) = sphericalbessely(m, x)
h1m(m, x) = sphericalbesselj(m, x) - im * sphericalbessely(m, x)
h2m(m, x) = sphericalbesselj(m, x) + im * sphericalbessely(m, x)
sph_jm(m, x) = x * jm(m, x)
sph_ym(m, x) = x * ym(m, x)
sph_h1m(m, x) = x * h1m(m, x)
sph_h2m(m, x) = x * h2m(m, x)
sph_jm_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_jm, l, x) 
sph_ym_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_ym, l, x) 
sph_h1m_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_h1m, l, x) 
sph_h2m_with_derivatives(l, x) = sph_bessel_with_derivatives(sph_h2m, l, x) 




h1m(m, x) = sphericalbesselj(m, x) - im * sphericalbessely(m, x)
h2m(m, x) = sphericalbesselj(m, x) + im * sphericalbessely(m, x)
@code_warntype h1m(3, 430.0)
@code_warntype h2m(3, 430.0)

@code_warntype sph_h1m(1, 430.0)
@code_warntype sph_h2m(1, 430.0)


typeof(sph_h2m_with_derivatives(1, 430.0))

@btime incident_wave($incident, $3, $(290.0), $(1.0))


k = wavenumber(f, μᵣ, εᵣ)

function incident_wave(incident, l, k, r) 
    if incident == true
        return sph_h1m_with_derivatives(l, k*r) 
    else
        return sph_h2m_with_derivatives(l, k*r)
    end
end

function fill_rs!(Rs, r_vec, lmax, k, incident)
    for l in 0:lmax
        for i in eachindex(r_vec) 
            x, y, z = r_vec[i]
            r = hypot(x, y, z)
            R, R´, R´´ = incident_wave(incident, l, k, r)
            Rs[i, l+1] = (R, R´, R´´)
        end
    end
end


basis = SphericalHarmonics(lmax)
Ylm, ∇Ylm = compute_with_gradients(basis, r_vec)

Ylm, ∇Ylm = AnalyticEMModes.testA(r_vec, lmax)

begin
    for l in 0:lmax
        for m in -l:l
            #idx = SpheriCart.lm2idx(l, m) # not public
            idx = l^2 + l + m + 1
            for ijk in eachindex(r_vec)
                A[ijk, idx] = te_sph_fields(r_vec[ijk], Rs[ijk, l+1], Ylm[ijk, idx], ∇Ylm[ijk, idx], k, μᵣ, εᵣ)
            end
        end
    end
end

return A



begin
    mode1 = myfields[:,13]
    fig = Figure()
    E = map(mode1, r_vec) do mi, ri
        E_r, E_θ, E_ϕ, H_r, H_θ, H_ϕ = mi
        r, θ, ϕ = cart2sph(ri...)
        spherical_to_cartesian_vector(0.0, H_θ, H_ϕ, θ, ϕ)
    end
    
    Ex = map(x->real(x[1]), E)
    Ey = map(x->real(x[2]), E)
    Ez = map(x->real(x[3]), E)
    
    Ef = (map((x,y,z) -> hypot(x,y,z), Ex, Ey, Ez))
    Ex ./= Ef
    Ey ./= Ef
    Ez ./= Ef
    factor = 0.15
    ax1 = Axis3(fig[1,1])
    arrows3d!(ax1, xcoords, ycoords, zcoords, Ex*factor, Ey*factor, Ez*factor, color = Ef, colormap = :jet)
    hidedecorations!(ax1)

    ax2 = Axis3(fig[1,2])
    Hρ = map(x->real(x[4]), mode1)
    mesh!(ax2, coords, conn, color = abs.(Hρ), colormap = :jet)
    hidedecorations!(ax2)
    fig
end



lmax =  5
basis = SphericalHarmonics(lmax)
Ylm, ∇Ylm = compute_with_gradients(basis, r_vec)



Rs = Matrix{NTuple{3, ComplexF64}}(undef, length(r_vec), lmax+1)

f, μᵣ, εᵣ = 10e9, 1.0, 1.0
k = wavenumber(f, μᵣ, εᵣ)
for l in 0:lmax
    for i in eachindex(r_vec) 
        x, y, z = r_vec[i]
        r = hypot(x, y, z)
        Rs[i, l+1] =  sph_h2m_with_derivatives(l, k*r)
    end
end


function spherical_to_cartesian_vector(gr, gθ, gϕ, θ, ϕ)
    # 1. Matriz de transformación T (o su transpuesta)
    # T[1,:] = (sinθcosϕ, cosθcosϕ, -sinϕ) 
    # T[2,:] = (sinθsinϕ, cosθsinϕ, cosϕ)
    # T[3,:] = (cosθ, -sinθ, 0)
    
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)

    # Componente X (gx): g_r * (r̂⋅x̂) + g_θ * (θ̂⋅x̂) + g_ϕ * (ϕ̂⋅x̂)
    gx = gr * (sinθ * cosϕ) + 
         gθ * (cosθ * cosϕ) - 
         gϕ * (sinϕ)

    # Componente Y (gy): g_r * (r̂⋅ŷ) + g_θ * (θ̂⋅ŷ) + g_ϕ * (ϕ̂⋅ŷ)
    gy = gr * (sinθ * sinϕ) + 
         gθ * (cosθ * sinϕ) + 
         gϕ * (cosϕ)

    # Componente Z (gz): g_r * (r̂⋅ẑ) + g_θ * (θ̂⋅ẑ) + g_ϕ * (ϕ̂⋅ẑ)
    gz = gr * (cosθ) - 
         gθ * (sinθ) 
         # + gϕ * 0
         
    return (gx, gy, gz)
end