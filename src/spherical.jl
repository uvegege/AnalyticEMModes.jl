#=
ψₘₙ(r, θ, ϕ) = R(r)Θ(θ)Φ(ϕ)

d/dr(r^2dR/dr) + [k^2*r^2 - n(n+1)]R = 0 -> Spherical Bessel Equation
1/sin(θ) d/dθ(sinθ * dΘ/dθ) + [n(n+1) - m^2/sin^2(θ)]Θ = 0 -> Legendre Equation

ψₘₙ(r, θ, ϕ) = [aₙjₙ(kcr) + bₙyₙ(kcr)][cₘₙPₙᵐ(cosθ) + dₘₙQₙᵐ(cosθ)][eₘcos(mϕ) + fₘsin(mϕ)]

=#
