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