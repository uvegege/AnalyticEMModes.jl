# Cylindrical Coordinate System

## Metric and Coordinate Description

The cylindrical coordinate system $(r, \theta, z)$ is characterized by the metric:

```math
ds^2 = dr^2 + r^2 d\theta^2 + dz^2
```

with scale factors:
```math
h_r = 1, \quad h_\theta = r, \quad h_z = 1
```

## Scalar Helmholtz Equation

The scalar Helmholtz equation in cylindrical coordinates:

```math
\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial \psi}{\partial r}\right) + \frac{1}{r^2}\frac{\partial^2 \psi}{\partial \theta^2} + \frac{\partial^2 \psi}{\partial z^2} + k^2 \psi = 0
```

## Separation of Variables and Bessel Functions

Assuming a solution of the form $\psi(r, \theta, z) = R(r)\Theta(\theta)Z(z)$, we obtain:

```math
\frac{1}{rR}\frac{d}{dr}\left(r\frac{dR}{dr}\right) + \frac{1}{r^2\Theta}\frac{d^2\Theta}{d\theta^2} + \frac{1}{Z}\frac{d^2Z}{dz^2} + k^2 = 0
```

This leads to the separation equations:

**Angular equation:**
```math
\frac{d^2\Theta}{d\theta^2} + m^2 \Theta = 0 \quad \Rightarrow \quad \Theta(\theta) = e^{\pm j m\theta}
```
where $m = 0, 1, 2, \dots$ for single-valuedness.

**Radial equation (Bessel's equation):**
```math
r^2 \frac{d^2R}{dr^2} + r \frac{dR}{dr} + (k_c^2 r^2 - m^2)R = 0
```
where $k_c^2 = k^2 + \gamma^2$.

The solutions are Bessel functions:
```math
R(r) = A J_m(k_c r) + B Y_m(k_c r)
```

- $J_m(k_c r)$: Bessel function of first kind (finite at $r=0$)
- $Y_m(k_c r)$: Bessel function of second kind (singular at $r=0$)


## Potential Functions and Field Relations

### TM Modes (E-modes) - Electric Vector Potential

Let $\mathbf{F} = \mathbf{\hat{z}} F_z(r, \theta)e^{-\gamma z}$, where $F_z$ satisfies:
```math
\nabla_t^2 F_z + k_c^2 F_z = 0
```

The field components are:
```math
E_r = -\frac{1}{\epsilon} \frac{1}{r} \frac{\partial F_z}{\partial \theta}, \quad
E_\theta = \frac{1}{\epsilon} \frac{\partial F_z}{\partial r}, \quad
E_z = 0
```
```math
H_r = -\frac{\gamma}{j\omega\mu\epsilon} \frac{\partial F_z}{\partial r}, \quad
H_\theta = -\frac{\gamma}{j\omega\mu\epsilon} \frac{1}{r} \frac{\partial F_z}{\partial \theta}, \quad
H_z = \frac{k_c^2}{j\omega\mu\epsilon} F_z
```