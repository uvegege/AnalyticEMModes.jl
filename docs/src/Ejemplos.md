# Examples and Gallery

To facilitate visualization of analytical results and validate their accuracy, we compare them with numerical solutions obtained using the finite element method. For this purpose, we use **Gridap.jl** to solve the electromagnetic eigenvalue problem.

> **Note**: Numerical results obtained with Gridap are not executed during documentation generation for simplicity, but the code used to generate the images is included as reference.

## Dependencies for Numerical Comparison

To reproduce the numerical comparisons, the following packages are required:

```julia
using Gmsh
using Gmsh: gmsh
using Gridap
using GridapGmsh
using Gridap.Helpers
using Gridap.Arrays
using Gridap.CellData
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Visualization
using Arpack
using GLMakie
```

## Helper Functions

### Mesh Data Extraction

The following function extracts node coordinates and element connectivity from Gmsh-generated meshes:

```julia
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
```

### Eigenvalue Problem in Gridap

The electromagnetic eigenvalue problem can be formulated concisely in Gridap. This function sets up the weak form of the Helmholtz equation for either TE or TM modes:

```julia
function eigenmode_gridap(model, T; order = 1)
    reffe = ReferenceFE(lagrangian, Float64, order)
    if T == :TM
        V = FESpace(model, reffe; conformity=:H1, dirichlet_tags=["PEC"])
    else

        V = FESpace(model, reffe; conformity=:H1)
    end
    U = V
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2*order)

    a(u,v) = ∫(∇(u)⋅∇(v))dΩ
    b(u,v) = ∫(u*v)dΩ

    A = assemble_matrix(a,U,V)
    B = assemble_matrix(b,U,V)
    return A, B, U, V
end
```

The eigenvalues λ obtained from this formulation are related to the cutoff wavenumber by kc² = λ. For TE modes, Neumann boundary conditions are applied (natural boundary condition in H1 space), while TM modes require Dirichlet boundary conditions (Ez = 0 at PEC walls).

## Examples by Geometry

Detailed examples for each geometry are organized in the following sections:

* **[Circular Waveguides](Cylindrical/Circular.md)**: Comparison of TE and TM modes in circular waveguides
* **[Coaxial Waveguides](Cylindrical/Coaxial.md)**: Visualization of TE, TM, and TEM modes in coaxial structures
* **[Radial and Wedge Waveguides](Cylindrical/Radial.md)**: Modes in radial waveguides and wedge structures
* **[Elliptic Waveguides](Elliptic/elliptic.md)**: Modes using Mathieu functions
* **[Rectangular Waveguides](Rectangular/rectangular.md)**: Analytical solutions in Cartesian coordinates
