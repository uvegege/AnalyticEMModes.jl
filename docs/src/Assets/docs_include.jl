#cd("C:\\MisProyecto\\Upload\\AnalyticEMModes\\docs")
#using Pkg
#Pkg.activate(;temp = true)
#Pkg.add("Gmsh")
#Pkg.add("Gridap")
#Pkg.add("Arpack")
#Pkg.add("GridapGmsh")

#using Gridap
#using GridapGmsh
#
#using Gridap.Helpers
#using Gridap.Arrays
#using Gridap.CellData
#using Gridap.Geometry
#using Gridap.ReferenceFEs
#using Gridap.Visualization
#
#using Arpack
using Gmsh
using Gmsh: gmsh
using GLMakie
include(raw"C:\MisProyecto\Upload\AnalyticEMModes\src\AnalyticEMModes.jl")
#using .AnalyticEMModes

include(joinpath(@__DIR__, "example_mesh.jl"))

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
    Connectivity = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    nodeTags, nodeCoords, _ = Gmsh.gmsh.model.mesh.getNodes()

    coordinates = reshape(nodeCoords, (3, length(nodeTags)))
    Gmsh.gmsh.finalize()
    return coordinates, Connectivity'
end

function mesh_data_radial(name)
    Gmsh.gmsh.initialize()
    Gmsh.gmsh.open(name)
    Gmsh.gmsh.model.mesh.renumberNodes()
    Gmsh.gmsh.model.mesh.renumberElements()

    dt = Gmsh.gmsh.model.getPhysicalGroups(2)[1]
    t = Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dt...)[1]
    element_tags, nodeTagsPerElement = gmsh.model.mesh.getElementsByType(2, t)

    tri_connectivity_flat = nodeTagsPerElement
    num_nodes_per_element = 3 # Triangles
    num_elements = length(tri_connectivity_flat) ÷ num_nodes_per_element
    connectivity = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    nodeTags, nodeCoords, _ = Gmsh.gmsh.model.mesh.getNodes()
    coordinates = reshape(nodeCoords, (3, length(nodeTags)))
    
    dt = Gmsh.gmsh.model.getPhysicalGroups(2)[2]
    t = reduce(vcat, Gmsh.gmsh.model.getEntitiesForPhysicalGroup(dt...))
    et = gmsh.model.mesh.getElementsByType.(2, t)
    tri_connectivity_flat = reduce(vcat, getindex.(et, 2))
    num_elements = length(tri_connectivity_flat) ÷ num_nodes_per_element
    connectivity2 = reshape(tri_connectivity_flat, (num_nodes_per_element, num_elements))
    Gmsh.gmsh.finalize()

    return coordinates, connectivity', connectivity2'
end

