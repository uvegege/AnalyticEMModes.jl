# Power Normalization

The normalization constant is derived from the power flow condition:

```math
P = \frac{1}{2} \text{Re} \int_0^a \int_0^b (\mathbf{E} \times \mathbf{H}^*) \cdot \mathbf{\hat{z}}  dx dy = 1 \text{ Watt}
```

## Power normalization in Rectangular Waveguides

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


## Power normalization in Circular and Coaxial Waveguides

### Field Representation

Let's define the potential function as 

```math
\psi(\rho) = A_m J_m(k_c\rho) + B_m Y_m(k_c\rho)
```
with derivative

```math
\frac{\partial\psi}{\partial \rho} = \psi^{\prime}(\rho) = k_c \left( A_m J_m'(k_c\rho) + B_m Y_m'(k_c\rho) \right)
```

For **TM modes**, the longitudinal electric field is:

```math
H_z = \psi(\rho) \cos(m\theta)
```

The transverse field components are:

```math
E_\rho = -\frac{j\beta}{k_c^2} \psi^{\prime}(\rho) \cos(m\theta), \quad  E_\theta = \frac{j\beta}{k_c^2} \frac{m}{\rho} \psi(\rho) \sin(m\theta)
```

```math
H_\rho = -\frac{j\omega\varepsilon}{k_c^2} \frac{m}{\rho} \psi(\rho) \sin(m\theta), \quad H_\theta = -\frac{j\omega\varepsilon}{k_c^2} \psi^{\prime}(\rho) \cos(m\theta)
```

### Power Flow 

The time-average power density is:

```math
S_z = \frac{1}{2} \left( E_\rho H_\theta^* - E_\theta H_\rho^*\right)
```

Substituting the field components:

```math
S_z = \frac{1}{2}\frac{\beta\omega\varepsilon}{k_c^4} \left( {\psi^\prime}^2\text{cos}^2\left( m\theta\right) + \frac{m^2}{\rho^2}{\psi}^2\text{sin}^2\left( m\theta\right) \right)
```

The total power is

```math
P = \frac{1}{2}\frac{\beta\omega\varepsilon}{k_c^4}\int_b^a \int_0^{2\pi} \left( {\psi^\prime}^2\text{cos}^2\left( m\theta\right) + \frac{m^2}{\rho^2}{\psi}^2\text{sin}^2\left( m\theta\right) \right) \,\rho\, d\rho d\theta
```

### Angular integration

The angular integrals depend on the azimuthal mode number m:

```math
\int_0^{2\pi} \cos^2(m\theta) \, d\theta = \begin{cases}
\pi & \text{if } m \geq 1 \\
2\pi & \text{if } m = 0
\end{cases}, \quad
\int_0^{2\pi} \sin^2(m\theta) \, d\theta = \begin{cases}
\pi & \text{if } m \geq 1 \\
0 & \text{if } m = 0
\end{cases}
```

Thus, the power becomes:

```math
P = \frac{1}{2} \frac{\beta\omega\varepsilon}{k_c^4} \cdot \pi \cdot I
```

where 

```math
I = \int_b^a \left[ (\psi^{\prime})^2 + \frac{m^2}{\rho^2} \psi^2 \right] \rho \, d\rho
```

### Simplifying the Radial Integral

The integral **I** can be simplified using multiple identities and relations. Consider:

```math
F = \psi
```

```math
\frac{\partial}{\partial\rho}\left( \rho F F^\prime  \right) = FF^\prime + \rho F^\prime F^\prime + \rho F F^{\prime\prime} = FF^\prime + \rho \left(F^\prime\right)^2+ \rho F F^{\prime\prime}
```

Rewriting the above expression:

```math
\rho \left(F^\prime\right)^2 = \frac{\partial}{\partial\rho}\left( \rho F F^\prime  \right) - \rho F F^{\prime\prime} - F F^\prime 
```

Integrating on both sides:

```math
\int_b^a \rho \left(F^\prime\right)^2  d\rho = \left[\rho F F^\prime\right]_b^a - \int_b^a F F^\prime d\rho - \int_b^a\rho F F^{\prime\prime} d\rho
```

We want to simplify the term with the second derivative. If we multiply Bessel's ODE we get:

```math
\left( \rho^2 F^{\prime\prime} + \rho F^\prime + \left( k_c^2\rho^2 - m^2\right)F\right)\rho F = \left(0\right) \rho F
```
```math
\rho F F^{\prime\prime} + F F^{\prime} + \left( k_c^2 \rho - \frac{m^2}{\rho}\right)F^2 = 0
```
```math
\rho F F^{\prime\prime} = -F F^{\prime} - \left(k_c^2\rho-\frac{m^2}{\rho}\right)F^2
```

After algebraic manipulation with the previous results, we get:

```math
\int_b^a \rho \left(F^\prime\right)^2  d\rho = \left[\rho F F^\prime\right]_b^a - \int_b^a F F^\prime d\rho - \int_b^a \left( -F F^{\prime} - \left(k_c^2\rho-\frac{m^2}{\rho}\right)F^2\right) d\rho
```

```math
\int_b^a \rho \left(F^\prime\right)^2  d\rho = \left[\rho F F^\prime\right]_b^a - \int_b^a F F^\prime d\rho + \int_b^a F F^{\prime}d\rho +\int_b^a  \left(k_c^2\rho-\frac{m^2}{\rho}\right)F^2 d\rho
```

```math
\int_b^a \rho \left(F^\prime\right)^2  d\rho = \left[\rho F F^\prime\right]_b^a +\int_b^a k_c^2\rho F^2 d\rho - \int_b^a \frac{m^2}{\rho}F^2 d\rho
```

If we reagroup the terms we can get exactly $I$ on the right hand side

```math
I = \int_b^a \rho \left(F^\prime\right)^2  d\rho + \int_b^a \frac{m^2}{\rho}F^2 d\rho = \left[\rho F F^\prime\right]_b^a + k_c^2\int_b^a \rho F^2 d\rho
```

Now we just need to solve the last part of the rhs:

```math
k_c^2\int_b^a \rho F^2 d\rho
```

### Bessel Function Integrals

The function $F$ depends of $\rho$. With the change of variable $x = k_c \rho$ we get

```math
\frac{k_c^2}{k_c^2}\int_{x_b}^{x_a} x F(x)^2 dx
```

The remaining integral can be evaluated using the identity:

```math
\int z C_\mu(az)D_\mu(az) = \frac{1}{4} z^2 \left(  2 C_\mu(az)D_\mu(az) - C_{\mu-1}(az)D_{\mu+1}(az) - C_{\mu+1}(az)D_{\mu-1}(az)\right)
```

where ``C_\mu`` and ``D_\mu`` are cylindrical functions(``J_\mu``, ``Y_\mu`` or ``H_\mu``). 

For ``m = 0`` the result is almost identical. The derivation follows the same steps, yielding:

```math
I = \int_b^a \rho \left(F^\prime\right)^2  d\rho = \left[\rho F F^\prime\right]_b^a + k_c^2\int_b^a \rho F^2 d\rho
```

which is exactly the same radial integral result. The only difference appears in the angular integration: instead of multiplying by ``\pi``, we now multiply by ``2\pi``.


## Circular Waveguide

For circular waveguides (``b = 0``, ``B_m = 0``), the boundary term vanishes at ``\rho = 0`` and 

```math
\psi = Am J_m(k_c\rho)
```

```math
I = k_c a J_m(k_c a) J_m'(k_c a) + \frac{(k_c a)^2}{2} \left[ J_m^2(k_c a) - J_{m+1}(k_c a) J_{m-1}(k_c a) \right]
```

### TM Modes

```math
F_0 = \sqrt{ \frac{2}{\displaystyle \frac{\omega\mu\beta}{k_c^4} \pi I } }
```

When ``m=0`` 

```math
F_0 = \sqrt{ \frac{2}{\displaystyle \frac{\omega\mu\beta}{k_c^4} 2\pi I } }
```

### TE Modes

```math
F_0 = \sqrt{ \frac{2}{\displaystyle \frac{\omega\varepsilon\beta}{k_c^4} \pi I } }
```

When ``m=0`` 

```math
F_0 = \sqrt{ \frac{2}{\displaystyle \frac{\omega\varepsilon\beta}{k_c^4} 2\pi I } }
```

## Coaxial Waveguide

For coaxial waveguides, both Bessel functions are required since the domain excludes the origin:

### TM Modes
```math
B_m = -\frac{J_m(k_c b)}{Y_m(k_c b)}
```
### TE Modes
```math
B_m = -\frac{J_m'(k_c b)}{Y_m'(k_c b)}
```

The radial function becomes
```math
F = A_m J_m(k_c\rho) + B_m Y_m(k_c\rho)
```
```math
F^{\prime} = k_c \left[ A_m J_m'(k_c\rho) + B_m Y_m'(k_c\rho) \right]
```

The integral **I** has two parts:

```math
I_1 = \left[ xF(x)F^{\prime}(x)\right]_b^a = \left( a F(a) F^{\prime}(a)\right) - \left( b F(b) F^{\prime}(b)\right) 
```

and

```math
I_2 =\int_b^a x F^2(x) dx = \int_b^a A_m^2J_m^2(x) + B_m^2 Y_m^2(x) + 2A_mB_mJ_m(x)Y_m(x) dx
```

```math
I_2 =\int_b^a x F^2(x) dx = A_m^2\int_b^a J_m^2(x)dx + B_m^2 \int_b^aY_m^2(x)dx + 2A_mB_m \int_b^a J_m(x)Y_m(x) dx
```

All this integrals can be solved with closed form solutions.

The normalization constant ``F_0`` has exactly the same expression as in circular waveguides, only changing the computation of **I**.

## Radial and Wedge Waveguide

In radial and wedge waveguides, power propagates radially along the $\rho$ direction.

### Power Flow Components

The radial component of the Poynting vector is:

```math
S_\rho = \frac{1}{2} \left( E_\phi H_z^* - E_z H_\phi^*\right)
```

* **TE modes** ``E_z = 0`` so ``S_\rho = \frac{1}{2} E_\phi H_z^*``

* **TM modes** ``H_z = 0`` so ``S_\rho = \frac{1}{2} E_z H_\phi^*``.

### Field Representation

Let's define the field functions:

```math
\psi(\rho) = A_{mn}H_m^{(2)}(k_c\rho)
```
```math
\phi_m(\varphi) = \begin{cases} \cos(m\varphi) \\ \sin(m\varphi) \end{cases}
```
```math
Z_n(z)  = \begin{cases} \sin(\frac{n\pi}{h}z) \\ \cos(\frac{n\pi}{h}z) \end{cases}
```

where the longitudinal field ``E_z`` or ``H_z`` (*TM* or *TE* modes) is defined as ``\psi(\rho)\phi_m(\varphi)Z_n(z) ``.

## Power Integration

For **TE modes**, integrating the power flow:

```math
P = \frac{1}{2} \int_0^{2\pi} \int_0^h \frac{j\mu\omega}{k_c^2}\psi^{\prime}\text{cos}(m\varphi)\text{sin}(\frac{n\pi}{h}z) \psi^*\text{sin}(\frac{n\pi}{h}z)\text{cos}(m\varphi) \rho d\varphi dz
```

The longitudinal and angular integrals yield:

```math
\int_0^h \text{sin}(\frac{n\pi}{h}z)  dz = \frac{h}{2}
```

```math
\int_0^{2\pi}\text{cos}^2(m\varphi) d\varphi = \alpha_m =\begin{cases} \pi & \text{if } m \geq 1 \\ 2\pi & \text{if } m = 0 \end{cases}
```

Thus, the total power simplifies to:

```math
P = \frac{1}{2} \text{Re} \left[ \frac{j\mu\omega}{k_c^2} \cdot \frac{h}{2} \cdot \alpha_m \cdot \rho \cdot \psi(\rho) \psi^{\prime}(\rho) \right]
```

```math
F_0 = \sqrt\frac{1}{P}
```

The same procedure and result applies to TM modes, replacing ``\mu`` by ``\varepsilon``.



