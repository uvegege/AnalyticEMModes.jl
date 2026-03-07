using LinearAlgebra
using StaticArrays
using Gmsh


function mesh_data(name)
    Gmsh.gmsh.initialize()
    Gmsh.gmsh.open(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()

    dt = Gmsh.gmsh.model.getPhysicalGroups(2)[1]
    t = Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dt...)[1]
    _, nodeTagsPerElement = gmsh.model.mesh.getElementsByType(2, t)

    tri_connectivity_flat = nodeTagsPerElement
    num_nodes_per_element = 3 # Triangles
    num_elements = length(tri_connectivity_flat) ÷ num_nodes_per_element
    connectivity = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    nodeTags, nodeCoords, _ = Gmsh.gmsh.model.mesh.getNodes()

    coordinates = reshape(nodeCoords, (3, length(nodeTags)))
    Gmsh.gmsh.finalize()
    return coordinates, connectivity'
end

# ---- utilities: read vertex indices for element e robustly ----
@inline function _tri_nodes(conn, e)
    # conn can be:
    #  - Matrix{Int}: (3, Ne) or (Ne, 3)
    #  - Vector{NTuple{3,Int}} or Vector{SVector{3,Int}} or Vector{Vector{Int}}
    if conn isa AbstractMatrix
        if size(conn,1) == 3
            return (conn[1,e], conn[2,e], conn[3,e])
        elseif size(conn,2) == 3
            return (conn[e,1], conn[e,2], conn[e,3])
        else
            error("conn matrix must be 3×Ne or Ne×3")
        end
    else
        ce = conn[e]
        return (ce[1], ce[2], ce[3])
    end
end

@inline function _pt(coord, i)
    # coord is dim×N (dim=2 or 3)
    if size(coord,1) == 2
        return SVector(coord[1,i], coord[2,i], 0.0)
    elseif size(coord,1) == 3
        return SVector(coord[1,i], coord[2,i], coord[3,i])
    else
        error("coord must be 2×N or 3×N")
    end
end

@inline function tri_area(p1, p2, p3)
    return 0.5 * norm(cross(p2-p1, p3-p1))
end

"""
    integrate_mesh_tri(coord, conn; f, quad=:tri1)

Integrate scalar function f(x::SVector{3}) over a triangular surface mesh.
Returns ∑_tri ∫_tri f dS.

quad:
- :tri1  -> 1-point (centroid), weight = area
- :tri3  -> 3-point degree-2 (barycentric), weights = area/3
"""
function integrate_mesh_tri(coord, conn; f::F, quad::Symbol=:tri1) where F
    Ne = conn isa AbstractMatrix ? (size(conn,1)==3 ? size(conn,2) : size(conn,1)) : length(conn)
    acc = zero(f(_pt(coord,1)))  # works for Float/Complex
    for e in 1:Ne
        i1,i2,i3 = _tri_nodes(conn, e)
        p1 = _pt(coord, i1); p2 = _pt(coord, i2); p3 = _pt(coord, i3)
        A = tri_area(p1,p2,p3)
        if A == 0
            continue
        end
        if quad === :tri1
            xc = (p1+p2+p3)/3
            acc += f(xc) * A
        elseif quad === :tri3
            # barycentric 3-pt rule
            # points: (2/3,1/6,1/6) permutations, weights = A/3
            x1 = (2p1 + p2 + p3)/4
            x2 = (p1 + 2p2 + p3)/4
            x3 = (p1 + p2 + 2p3)/4
            w = A/3
            acc += w*(f(x1) + f(x2) + f(x3))
        else
            error("quad must be :tri1 or :tri3")
        end
    end
    return acc
end