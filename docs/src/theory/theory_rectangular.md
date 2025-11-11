# Cartesian Coordinate System

## Metric and Coordinate Description

The Cartesian coordinate system $(x, y, z)$ is characterized by the metric:

```math
ds^2 = dx^2 + dy^2 + dz^2
```

with scale factors:
```math
h_x = 1, \quad h_y = 1, \quad h_z = 1
```

## Scalar Helmholtz Equation

The scalar Helmholtz equation in Cartesian coordinates:

```math
\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2} + \frac{\partial^2 \psi}{\partial z^2} + k^2 \psi = 0
```

## Separation of Variables

Assuming a solution of the form $\psi(x, y, z) = X(x)Y(y)Z(z)$, we obtain:

```math
\frac{1}{X}\frac{d^2X}{dx^2} + \frac{1}{Y}\frac{d^2Y}{dy^2} + \frac{1}{Z}\frac{d^2Z}{dz^2} + k^2 = 0
```

This leads to the separation equations:
```math
\frac{d^2X}{dx^2} + k_x^2 X = 0, \quad
\frac{d^2Y}{dy^2} + k_y^2 Y = 0, \quad
\frac{d^2Z}{dz^2} + k_z^2 Z = 0
```
with the constraint:
```math
k_x^2 + k_y^2 + k_z^2 = k^2
```

## Potential Functions and Field Solutions

We define vector potentials:

### TM Modes - Electric Vector Potential

Let ($\mathbf{F} = \mathbf{\hat{z}} F_z(x, y)e^{-\gamma z}$), where ($F_z$) satisfies:
```math
\nabla_t^2 F_z + k_c^2 F_z = 0
```
with ($k_c^2 = k_x^2 + k_y^2$).

The field components are:
```math
E_x = -\frac{1}{\epsilon} \frac{\partial F_z}{\partial y}, \quad
E_y = \frac{1}{\epsilon} \frac{\partial F_z}{\partial x}, \quad
E_z = 0
```
```math
H_x = -\frac{1}{j\omega\mu\epsilon} \frac{\partial^2 F_z}{\partial x \partial z}, \quad
H_y = -\frac{1}{j\omega\mu\epsilon} \frac{\partial^2 F_z}{\partial y \partial z}, \quad
H_z = \frac{1}{j\omega\mu\epsilon} (k^2 + \frac{\partial^2}{\partial z^2}) F_z
```

### TE Modes - Magnetic Vector Potential  

Let ($\mathbf{A} = \mathbf{\hat{z}} A_z(x, y)e^{-\gamma z}$), where ($A_z$) satisfies:
```math
\nabla_t^2 A_z + k_c^2 A_z = 0
```

The field components are:
```math
E_x = -\frac{1}{j\omega\mu\epsilon} \frac{\partial^2 A_z}{\partial x \partial z}, \quad
E_y = -\frac{1}{j\omega\mu\epsilon} \frac{\partial^2 A_z}{\partial y \partial z}, \quad
E_z = \frac{1}{j\omega\mu\epsilon} (k^2 + \frac{\partial^2}{\partial z^2}) A_z
```
```math
H_x = \frac{1}{\mu} \frac{\partial A_z}{\partial y}, \quad
H_y = -\frac{1}{\mu} \frac{\partial A_z}{\partial x}, \quad
H_z = 0
```

# Rectangular Waveguide Solutions

For a rectangular waveguide with PEC walls at \(x = 0, a\) and \(y = 0, b\), the general solution to the Helmholtz equation is:

```math
\psi(x, y) = [A\cos(k_x x) + B\sin(k_x x)] \cdot [C\cos(k_y y) + D\sin(k_y y)]
```

Applying the PEC boundary conditions:

- **TM modes**: Require ($\psi = 0$) on all boundaries, leading to sine functions in both directions
- **TE modes**: Require normal derivatives vanishing on boundaries, leading to cosine functions

This yields the well-known modal solutions:


### TM Modes

```math
F_z(x, y) = \sin\left(\frac{m\pi x}{a}\right) \sin\left(\frac{n\pi y}{b}\right), \quad m, n = 1, 2, 3, \dots
```

Cutoff wavenumber:

```math
k_c = \sqrt{\left(\frac{m\pi}{a}\right)^2 + \left(\frac{n\pi}{b}\right)^2}
```

### TE Modes

```math
A_z(x, y) = \cos\left(\frac{m\pi x}{a}\right) \cos\left(\frac{n\pi y}{b}\right), \quad m, n = 0, 1, 2, \dots
```

but not $m = n = 0$.

Cutoff wavenumber:

```math
k_c = \sqrt{\left(\frac{m\pi}{a}\right)^2 + \left(\frac{n\pi}{b}\right)^2}
```


## Power Normalization

The normalization constant is derived from the power flow condition:

```math
P = \frac{1}{2} \text{Re} \int_0^a \int_0^b (\mathbf{E} \times \mathbf{H}^*) \cdot \mathbf{\hat{z}}  dx dy = 1 \text{ Watt}
```

### For TE Modes

Using the field components (with the factors implemented in the package):
```math
	E_x = F_0 \frac{\omega\mu}{k_c^2}\frac{mx}{a}\sin \left(\frac{mx}{a}\right) \cos\left(\frac{ny}{b}\right), \quad
	E_y = -F_0 \frac{\omega\mu}{k_c^2}\frac{ny}{b}{\epsilon} F_0 \cos\left(\frac{mx}{a}\right)\sin\left(\frac{ny}{b}\right)
```

```math
	H_x = -F_0 \frac{\beta}{k_c^2}\frac{ny}{b}{\epsilon} F_0 \cos\left(\frac{mx}{a}\right)\sin\left(\frac{ny}{b}\right), \quad
	H_y = F_0 \frac{\beta}{k_c^2}\frac{mx}{a}\sin \left(\frac{mx}{a}\right) \cos\left(\frac{ny}{b}\right)
```

The z-component of the Poynting vector is:
```math
S_z = \frac{1}{2} (E_x H_y^* - E_y H_x^*)
```

Substituting and simplifying:
```math
S_z = \frac{1}{2} \frac{\beta F_0^2 \omega \mu}{k_c^4} \left[k_y^2 \sin^2(k_y y)\cos^2(k_x x) + k_x^2 \cos^2(k_y y)\sin^2(k_x x)\right]
```
The power is

```math
P_z = \int_{0}^a\int_0^b \frac{1}{2} \frac{\beta F_0^2 \omega \mu}{k_c^4} \left[k_y^2 \sin^2(k_y y)\cos^2(k_x x) + k_x^2 \cos^2(k_y y)\sin^2(k_x x)\right]
```

Integrating over the cross-section and using the orthogonality of trigonometric functions:
```math
\int_0^a \sin^2(k_x x) dx = \frac{a}{2}, \quad \int_0^a \cos^2(k_x x) dx = \frac{a}{2}
```
```math
\int_0^b \cos^2(k_y y) dy = \frac{b}{2}, \quad \int_0^b \sin^2(k_y y) dy = \frac{b}{2}
```

Defining:
```math
\delta_m = \begin{cases} 1 & \text{if } n \geq 1 \\ 2 & \text{if } n = 0 \end{cases}, \quad
\delta_n = \begin{cases} 1 & \text{if } m \geq 1 \\ 2 & \text{if } m = 0 \end{cases}
```
We get:
```math
P_z = \frac{1}{2} \frac{\beta F_0^2 \omega \mu}{k_c^4} \left( \int_{0}^a\int_0^b k_y^2 \sin^2(k_y y)\cos^2(k_x x) +  \int_{0}^a\int_0^b k_x^2 \cos^2(k_y y)\sin^2(k_x x) \right)
```

```math
P_z = \frac{1}{2} \frac{\beta F_0^2 \omega \mu}{k_c^4} \left( \delta_m k_y^2\frac{ab}{4} +  \delta_n k_x^2 \frac{ab}{4}\right) = = \frac{1}{2} \frac{\beta F_0^2 \omega \mu}{k_c^4} \left( \delta_m \frac{n^2\pi^2}{b^2} \frac{ab}{4} + \delta_n \frac{m^2\pi^2}{a^2}\frac{ab}{4} \right)
```

```math
P_z = \frac{1}{2} \frac{\beta F_0^2 \omega \mu \pi^2}{4k_c^4} \left( \delta_m n^2 \frac{a}{b} +  \delta_n m^2\frac{b}{a} \right)
```

Setting $P = 1$ and solving for $F_0$:
```math
F_0 = \sqrt{\frac{8k_c^4}{\beta \omega \mu \pi^2 \left( \delta_m n^2 \frac{a}{b} +  \delta_n m^2\frac{b}{a} \right)}}
```

For **TM** modes:

```math
F_0 = \sqrt{\frac{8k_c^4}{\beta \omega \varepsilon \pi^2 \left( \delta_m n^2 \frac{a}{b} +  \delta_n m^2\frac{b}{a} \right)}}
```

## References

Harrington, R. F. (1961). *Time-Harmonic Electromagnetic Fields*, Chapter 3.  
Stratton, J. A. (1941). *Electromagnetic Theory*, Chapter 5.  
Balanis, C. A. (2012). *Advanced Engineering Electromagnetics*, Chapter 8.
