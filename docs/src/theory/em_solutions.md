# Electromagnetic Modes Theory

## 1. What are Electromagnetic Modes?

In electromagnetic theory, many physical problems admit solutions that form *distinct spatial field patterns* satisfying both **Maxwell's equations** and specific **boundary or radiation conditions**.

These patterns are known as **electromagnetic modes** — self-consistent field configurations that can exist and evolve independently.

Mathematically, modes emerge as **eigenfunctions** of the vector Helmholtz equation under given constraints. Each mode is characterized by:

- A unique **field distribution** in space
- A specific **frequency characteristic** (cutoff, resonant, or propagation)
- A defined **propagation behavior** (traveling, standing, or evanescent waves)
- **Orthogonality** properties allowing linear superposition without interference

## 2. Helmholtz Equation and Modal Solutions

In source-free, linear, homogeneous regions, Maxwell’s equations can be reduced to a **vector Helmholtz equation**:

```math
\nabla^2 \mathbf{E} + k^2 \mathbf{E} = 0, \quad
\nabla^2 \mathbf{H} + k^2 \mathbf{H} = 0
```

where

```math
k = \omega \sqrt{\mu \, \varepsilon}
```

is the **wavenumber** in the medium.

When the geometry allows separation of variables (e.g., in rectangular, circular, or other orthogonal coordinates),  
each Cartesian or curvilinear component satisfies a **scalar Helmholtz equation** of the form:

```math
\nabla_t^2 \psi + k_c^2 \psi = 0
```

where $\psi$ is a potential function (often $E_z$ or $H_z$),  
and $k_c$ is the **cutoff wavenumber** determined by the boundary conditions.  

The remaining transverse field components ($E_x$, $E_y$, $H_x$, $H_y$) are then derived from $\psi$ using Maxwell’s curl relations.


## 3. Classification of Modes

Depending on which longitudinal component is zero, three fundamental mode types exist:

| Mode Type | Condition | Typical Example |
|------------|------------|----------------|
| **TE (Transverse Electric)** | \( E_z = 0 \) | Waveguides |
| **TM (Transverse Magnetic)** | \( H_z = 0 \) | Waveguides |
| **TEM (Transverse Electromagnetic)** | \( E_z = H_z = 0 \) | Coaxial or parallel-plate line |

> TEM modes require two or more conductors (or open boundaries) to support a non-zero voltage and current pair.

These classifications reflect the physical symmetry of the problem and simplify the field representation, since only one scalar potential (either \( H_z \) or \( E_z \)) needs to be solved.

Hybrid modes are electromagnetic solutions where both longitudinal field components ($E_z$ and $H_z$) are non-zero, typically appearing in dielectric waveguides, inhomogeneous structures, etc. While important in many practical applications like optical fibers, hybrid modes are not the primary focus of this package, which concentrates on canonical structures with pure TE/TM analytical solutions.


## 4. Propagation and Cutoff

For a given mode, the propagation constant $\gamma$ determines the mode behavior:

```math
\gamma = \sqrt{k_c^2 - k^2} = \begin{cases}
j\beta = j\sqrt{k^2 - k_c^2} & \text{for } k > k_c \text{ (propagating)} \\
\alpha = \sqrt{k_c^2 - k^2} & \text{for } k < k_c \text{ (evanescent)}
\end{cases}
```

The **cutoff frequency** is defined as:

```math
f_c = \frac{k_c}{2 \pi \sqrt{\mu \, \varepsilon}}
```

Below cutoff, modes decay exponentially as $e^{-\alpha z}$ and carry no real power.


## 5. Power Flow and the Poynting Vector

Once the modal electric and magnetic fields are known,  
the instantaneous **power density** is given by the **Poynting vector**:

```math
\mathbf{S} = \mathbf{E} \times \mathbf{H}^*
```

The time-averaged power flow through a surface is:

```math
P = \frac{1}{2} \, \text{Re} \int_S (\mathbf{E} \times \mathbf{H}^*) \cdot \hat{n} \, dS
```

In propagating modes, $P$ is positive and represents transmitted energy.  
In evanescent modes, ($P = 0$): the fields store energy locally but do not transport it.

## 6. Orthogonality and Normalization

Modal fields satisfy an orthogonality condition that makes them particulary useful in analytical and numerical analysis:

```math
\int_S (\mathbf{E}_m \times \mathbf{H}_n^*) \cdot \hat{\mathbf{n}}\, d\mathbf{S} = 0 \quad \text{for} \quad m \neq n```

This means each mode can propagate independently, and total fields can be expanded as a linear combination of modes:

```math
\mathbf{E} = \sum_m a_m \mathbf{E}_m, \quad \mathbf{H} = \sum_m a_m \mathbf{H}_m
```

## References

