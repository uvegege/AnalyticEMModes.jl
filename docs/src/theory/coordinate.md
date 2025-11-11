# Introduction to Orthogonal Coordinate Systems

## Mathematical Foundation

Analytic electromagnetic solutions rely on orthogonal coordinate systems where the vector Helmholtz equation is separable. In such systems, the wave equation admits solutions through the method of separation of variables.

### General Formulation

For any orthogonal coordinate system $(u_1, u_2, u_3)$ with scale factors $(h_1, h_2, h_3)$, the differential operators take the form:

#### Gradient

```math
\nabla \psi = \frac{1}{h_1} \frac{\partial \psi}{\partial u_1} \mathbf{\hat{u}_1} + 
\frac{1}{h_2} \frac{\partial \psi}{\partial u_2} \mathbf{\hat{u}_2} + 
\frac{1}{h_3} \frac{\partial \psi}{\partial u_3} \mathbf{\hat{u}_3}
```

#### Divergence

```math
\nabla \cdot \mathbf{A} = \frac{1}{h_1 h_2 h_3} \left[
\frac{\partial}{\partial u_1} (h_2 h_3 A_1) +
\frac{\partial}{\partial u_2} (h_3 h_1 A_2) +
\frac{\partial}{\partial u_3} (h_1 h_2 A_3)
\right]
```

#### Curl

```math
\nabla \times \mathbf{A} = \frac{1}{h_1 h_2 h_3} \begin{vmatrix}
h_1 \mathbf{\hat{u}_1} & h_2 \mathbf{\hat{u}_2} & h_3 \mathbf{\hat{u}_3} \\
\frac{\partial}{\partial u_1} & \frac{\partial}{\partial u_2} & \frac{\partial}{\partial u_3} \\
h_1 A_1 & h_2 A_2 & h_3 A_3
\end{vmatrix}
```


#### Laplacian


```math
\nabla^2 \psi = \frac{1}{h_1 h_2 h_3} \left[ 
\frac{\partial}{\partial u_1} \left( \frac{h_2 h_3}{h_1} \frac{\partial \psi}{\partial u_1} \right) +
\frac{\partial}{\partial u_2} \left( \frac{h_3 h_1}{h_2} \frac{\partial \psi}{\partial u_2} \right) +
\frac{\partial}{\partial u_3} \left( \frac{h_1 h_2}{h_3} \frac{\partial \psi}{\partial u_3} \right)
\right]
```

The vector Helmholtz equation for electromagnetic fields:
```math
\nabla^2 \mathbf{E} + k^2 \mathbf{E} = 0, \quad \nabla^2 \mathbf{H} + k^2 \mathbf{H} = 0
```

can be solved by assuming separable solutions of the form:

```math
\psi(u_1, u_2, u_3) = U_1(u_1) U_2(u_2) U_3(u_3)
```

## Coordinate Systems implemented

This package supports solutions in the following orthogonal coordinate systems:

* Cartesian: $(x, y, z)$ - for rectangular geometries
* Cylindrical: $(r, \theta, z)$ - for circular, coaxial, and wedge geometries
* Elliptical: $(\xi, \eta, z)$ - for elliptical waveguides
* Spherical: $(r, \theta, \phi)$ - for spherical radiators

Each coordinate system is particularly suited for specific geometries where boundary conditions align naturally with constant coordinate surfaces.

## Reference

Stratton, J. A. (1941). *Electromagnetic Theory*, Chapter 1. McGraw-Hill.

