# Spheroidal Wave Functions in Electromagnetic Theory
# Le-Wei Li, Xiao-Kang Kang, Mook-Seng Leong
# Appendix A


function Mˣ_oblate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (S * ∂R * sϕ * cmϕ - m * abs(ξ) / (factor1) * ∂S * R * cϕ * smϕ)
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (∂S * R * sϕ * cmϕ + m * η / factor3 * S * R * cϕ * smϕ)
    mϕ = 2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * cϕ * cmϕ

    return (mξ, mη, mϕ)
end

function Mˣ_oblate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (S * ∂R * sϕ * smϕ - m * abs(ξ) / (factor1) * ∂S * R * cϕ * (-cmϕ))
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (∂S * R * sϕ * smϕ + m * η / factor3 * S * R * cϕ * (-cmϕ))
    mϕ = 2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * cϕ * smϕ

    return (mξ, mη, mϕ)
end

function Mˣ_prolate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (S * ∂R * sϕ * cmϕ - m * abs(ξ) / (factor1) * ∂S * R * cϕ * smϕ)
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (∂S * R * sϕ * cmϕ + m * η / factor3 * S * R * cϕ * smϕ)
    mϕ = 2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * cϕ * cmϕ

    return (mξ, mη, mϕ)
end

function Mˣ_prolate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (S * ∂R * sϕ * smϕ - m * abs(ξ) / (factor1) * ∂S * R * cϕ * -cmϕ)
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (∂S * R * sϕ * smϕ + m * η / factor3 * S * R * cϕ * -cmϕ)
    mϕ = 2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * cϕ * smϕ

    return (mξ, mη, mϕ)
end

function Mʸ_oblate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (-S * ∂R * cϕ * cmϕ - m * abs(ξ) / (factor1) * S * R * sϕ * smϕ)
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (-∂S * R * cϕ * cmϕ + m * η / factor3 * S * R * sϕ * smϕ)
    mϕ = -2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * sϕ * cmϕ

    return (mξ, mη, mϕ)
end

function Mʸ_oblate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (-S * ∂R * cϕ * smϕ - m * abs(ξ) / (factor1) * S * R * sϕ * (-cmϕ))
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (-∂S * R * cϕ * smϕ + m * η / factor3 * S * R * sϕ * (-cmϕ))
    mϕ = -2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * sϕ * smϕ

    return (mξ, mη, mϕ)
end

function Mʸ_prolate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (-S * ∂R * cϕ * cmϕ - m * abs(ξ) / (factor1) * S * R * sϕ * smϕ)
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (-∂S * R * cϕ * cmϕ + m * η / factor3 * S * R * sϕ * smϕ)
    mϕ = -2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * sϕ * cmϕ

    return (mξ, mη, mϕ)
end

function Mʸ_prolate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = -2 * sqrt(factor1) / (d * sqrt(factor2)) * (-S * ∂R * cϕ * smϕ - m * abs(ξ) / (factor1) * S * R * sϕ * (-cmϕ))
    mξ = 2 * sqrt(factor3) / (d * sqrt(factor2)) * (-∂S * R * cϕ * smϕ + m * η / factor3 * S * R * sϕ * (-cmϕ))
    mϕ = -2 / (d * sqrt(factor2)) * (ξ * factor3 * ∂S * R + η * factor1 * S * ∂R) * sϕ * smϕ

    return (mξ, mη, mϕ)
end


function Mᶻ_oblate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = 2 * m * η / (d * sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mξ = 2 * m * ξ / (d * sqrt(factor2) * sqrt(factor1)) * S * R * (-smϕ)
    mϕ = 2 * sqrt(factor1) * sqrt(factor3) / (d * sqrt(factor2)) * (η * ∂S * R - abs(ξ) * S * ∂R) * cmϕ

    return (mξ, mη, mϕ)
end

function Mᶻ_oblate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = 2 * m * η / (d * sqrt(factor2) * sqrt(factor3)) * S * R * (-cmϕ)
    mξ = 2 * m * ξ / (d * sqrt(factor2) * sqrt(factor1)) * S * R * cmϕ
    mϕ = 2 * sqrt(factor1) * sqrt(factor3) / (d * sqrt(factor2)) * (η * ∂S * R - abs(ξ) * S * ∂R) * smϕ

    return (mξ, mη, mϕ)
end

function Mᶻ_prolate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = 2 * m * η / (d * sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mξ = 2 * m * ξ / (d * sqrt(factor2) * sqrt(factor1)) * S * R * (-smϕ)
    mϕ = 2 * sqrt(factor1) * sqrt(factor3) / (d * sqrt(factor2)) * (η * ∂S * R - abs(ξ) * S * ∂R) * cmϕ

    return (mξ, mη, mϕ)
end

function Mᶻ_prolate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = 2 * m * η / (d * sqrt(factor2) * sqrt(factor3)) * S * R * (-cmϕ)
    mξ = 2 * m * ξ / (d * sqrt(factor2) * sqrt(factor1)) * S * R * cmϕ
    mϕ = 2 * sqrt(factor1) * sqrt(factor3) / (d * sqrt(factor2)) * (η * ∂S * R - abs(ξ) * S * ∂R) * smϕ

    return (mξ, mη, mϕ)
end

function Mʳ_oblate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = m * abs(ξ) / (sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mξ = +m * η / (sqrt(factor2) * sqrt(factor1)) * S * R * smϕ
    mϕ = sqrt(factor1) * sqrt(factor3) / factor2 * (abs(ξ) * ∂S * R + η * S * ∂R) * cmϕ

    return (mξ, mη, mϕ)
end

function Mʳ_oblate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = m * abs(ξ) / (sqrt(factor2) * sqrt(factor3)) * S * R * (-cmϕ)
    mξ = +m * η / (sqrt(factor2) * sqrt(factor1)) * S * R * (-cmϕ)
    mϕ = sqrt(factor1) * sqrt(factor3) / factor2 * (abs(ξ) * ∂S * R + η * S * ∂R) * smϕ

    return (mξ, mη, mϕ)
end

function Mʳ_oblate_even(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = m * abs(ξ) / (sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mξ = -m * η / (sqrt(factor2) * sqrt(factor1)) * S * R * smϕ
    mϕ = sqrt(factor1) * sqrt(factor3) / factor2 * (abs(ξ) * ∂S * R - η * S * ∂R) * cmϕ

    return (mξ, mη, mϕ)
end

function Mʳ_oblate_odd(m, n, c, ξ, η, ϕ)

    d = 1
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    mη = m * abs(ξ) / (sqrt(factor2) * sqrt(factor3)) * S * R * (-cmϕ)
    mξ = -m * η / (sqrt(factor2) * sqrt(factor1)) * S * R * (-cmϕ)
    mϕ = sqrt(factor1) * sqrt(factor3) / factor2 * (abs(ξ) * ∂S * R - η * S * ∂R) * smϕ

    return (mξ, mη, mϕ)
end

#TODO: M⁺ and M⁻

# Nˣ functions

function Nˣ_oblate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (3 * η^2 + ξ^2 - 2) / factor2^2) # Oblate: (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((2 * ξ^2 + 1) - ξ^2) / (sqrt(factor1) * factor2^2) # Oblate: η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2)

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
             η * S * (v1 * ∂²R + ∂v1 * ∂R) -
             1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - m^2 * η * S * R / (factor3 * sqrt(factor1)) *
                                                                               cϕ * cmϕ + m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * sϕ * smϕ
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = -(η * factor3 * (η^2 + 3 * ξ^2 + 2)) / factor2^2 # Oblate: (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (1 - 2 * η^2) - η^2) / (sqrt(factor3) * factor2^2) # Oblate: (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2)

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R +
                                           1 / sqrt(factor3) * S * ∂R + factor1 * (
                                                                            v4 * ∂S + ∂v4 * S) * ∂R -
                                           m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * cϕ * cmϕ +
         m / (factor3) * S * (ξ / factor1 * R - ∂R) * sϕ * smϕ

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
        (1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R +
        1 / sqrt(factor3) * S * (
            sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * sϕ * cmϕ + m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S +
                                                 (-1) / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

function Nˣ_oblate_odd(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (3 * η^2 + ξ^2 - 2) / factor2^2) # Oblate: (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((2 * ξ^2 + 1) - ξ^2) / (sqrt(factor1) * factor2^2) # Oblate: η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2)

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
             η * S * (v1 * ∂²R + ∂v1 * ∂R) -
             1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - m^2 * η * S * R / (factor3 * sqrt(factor1)) *
                                                                               cϕ * smϕ + m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * sϕ * (-cmϕ)
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = -(η * factor3 * (η^2 + 3 * ξ^2 + 2)) / factor2^2 # Oblate: (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (1 - 2 * η^2) - η^2) / (sqrt(factor3) * factor2^2) # Oblate: (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2)

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R + 1 / sqrt(factor3) * S * ∂R + factor1 * (
                                                                                                            v4 * ∂S + ∂v4 * S) * ∂R - m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * cϕ * smϕ + m / (factor3) * S * (ξ / factor1 * R - ∂R) * sϕ * (-cmϕ)

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
        (1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R + 1 / sqrt(factor3) * S * (
        sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * sϕ * smϕ + m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S + 
        -1 / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

function Nˣ_prolate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2) # Oblate: (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2) # Oblate: η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2)

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
             η * S * (v1 * ∂²R + ∂v1 * ∂R) -
             1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - m^2 * η * S * R / (factor3 * sqrt(factor1)) *
                                                                               cϕ * cmϕ + m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * sϕ * smϕ
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2 # Oblate: (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2) # Oblate: (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2)

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R +
                                           1 / sqrt(factor3) * S * ∂R + factor1 * (
                                                                            v4 * ∂S + ∂v4 * S) * ∂R -
                                           m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * cϕ * cmϕ +
         m / (factor3) * S * (ξ / factor1 * R - ∂R) * sϕ * smϕ

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
        (1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R +
        1 / sqrt(factor3) * S * (
            sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * sϕ * cmϕ + m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S +
                                                 (-1) / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

function Nˣ_prolate_odd(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 =  (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2) # Oblate: (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2) # Oblate: η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2)

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
             η * S * (v1 * ∂²R + ∂v1 * ∂R) -
             1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - m^2 * η * S * R / (factor3 * sqrt(factor1)) *
                                                                               cϕ * smϕ + m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * sϕ * (-cmϕ)
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2 # Oblate: (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2) # Oblate: (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2)

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R + 1 / sqrt(factor3) * S * ∂R + factor1 * (
                                                                                                            v4 * ∂S + ∂v4 * S) * ∂R - m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * cϕ * smϕ + m / (factor3) * S * (ξ / factor1 * R - ∂R) * sϕ * (-cmϕ)

    #------
    #mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
    #    (1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R + 1 / sqrt(factor3) * S * (
    #    sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * sϕ * smϕ + m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S + 
    #    -1 / ((factor3)^(3 / 2)) * ∂S) * R -
    #    1 / sqrt(factor3) * S * (ξ / sqrt(factor1) * ∂R + 1 / (sqrt(factor1)^(3 / 2)) * R)) * cϕ * (-cmϕ)
    #)

    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
        (1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R + 1 / sqrt(factor3) * S * (
        sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * sϕ * smϕ + m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S + 
        -1 / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

# Nʸ functions

function Nʸ_oblate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (3 * η^2 + ξ^2 - 2) / factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((2 * ξ^2 + 1) - ξ^2) / (sqrt(factor1) * factor2^2) 

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
            η * S * (v1 * ∂²R + ∂v1 * ∂R) -
            1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - 
            m^2 * η * S * R / (factor3 * sqrt(factor1)) *
            sϕ * cmϕ - m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * cϕ * smϕ
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = -(η * factor3 * (η^2 + 3 * ξ^2 + 2)) / factor2^2 
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (1 - 2 * η^2) - η^2) / (sqrt(factor3) * factor2^2) 

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R +
                1 / sqrt(factor3) * S * ∂R + factor1 * (
                v4 * ∂S + ∂v4 * S) * ∂R -
                m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * sϕ * cmϕ -
                m / (factor3) * S * (ξ / factor1 * R - ∂R) * cϕ * smϕ

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
            -(1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R +
            1 / sqrt(factor3) * S * (
            sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * cϕ * cmϕ + 
            m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S +
            -1 / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

function Nʸ_oblate_odd(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (3 * η^2 + ξ^2 - 2) / factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((2 * ξ^2 + 1) - ξ^2) / (sqrt(factor1) * factor2^2) 

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
            η * S * (v1 * ∂²R + ∂v1 * ∂R) -
            1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - 
            m^2 * η * S * R / (factor3 * sqrt(factor1)) *
            sϕ * smϕ - m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * cϕ * (-cmϕ)
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = -(η * factor3 * (η^2 + 3 * ξ^2 + 2)) / factor2^2 
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (1 - 2 * η^2) - η^2) / (sqrt(factor3) * factor2^2) 

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R +
                1 / sqrt(factor3) * S * ∂R + factor1 * (
                v4 * ∂S + ∂v4 * S) * ∂R -
                m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * sϕ * smϕ -
                m / (factor3) * S * (ξ / factor1 * R - ∂R) * cϕ * (-cmϕ)

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
            -(1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R +
            1 / sqrt(factor3) * S * (
            sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * cϕ * smϕ + 
            m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S +
            -1 / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

function Nʸ_prolate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2) 
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2) 

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
            η * S * (v1 * ∂²R + ∂v1 * ∂R) -
            1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - 
            m^2 * η * S * R / (factor3 * sqrt(factor1)) * sϕ * cmϕ - 
            m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * cϕ * smϕ
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2 
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2)

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R +
            1 / sqrt(factor3) * S * ∂R + factor1 * (
            v4 * ∂S + ∂v4 * S) * ∂R -
            m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * sϕ * cmϕ -
            m / (factor3) * S * (ξ / factor1 * R - ∂R) * cϕ * smϕ

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
        -(1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R +
        1 / sqrt(factor3) * S * (
        sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * cϕ * cmϕ + 
        m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S +
        -1 / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

function Nʸ_prolate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    #------
    v1 = (factor1^(3 / 2) / factor2)
    ∂v1 = (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2) # Oblate: (ξ * sqrt(factor1) * (-3*η^2 + ξ^2 + 2) /factor2^2)
    v2 = ξ * sqrt(factor1) / factor2
    ∂v2 = η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2) # Oblate: η^2 * ((1-2*ξ^2) + ξ^2)/(sqrt(factor1) * factor2^2)

    mη = 4 / (k * d^2 * sqrt(factor2)) * (
            η * S * (v1 * ∂²R + ∂v1 * ∂R) -
            1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (∂v2 * R + v2 * ∂R) - 
            m^2 * η * S * R / (factor3 * sqrt(factor1)) * sϕ * smϕ - 
            m / sqrt(factor1) * (∂S + η / factor3 * S)) * R * cϕ * (-cmϕ)
    #------
    v3 = factor3^(3 / 2) / factor2
    ∂v3 = (η * factor3 * (η^2 - 3*ξ^2 + 2)) / factor2^2 
    v4 = η * sqrt(factor3) / factor2
    ∂v4 = (ξ^2 * (2*η^2 - 1) - η^2) / (sqrt(factor3) * factor2^2)

    mξ = -4 / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂²S + ∂v3 * ∂S) * R +
            1 / sqrt(factor3) * S * ∂R + factor1 * (
            v4 * ∂S + ∂v4 * S) * ∂R -
            m^2 * ξ * S * R / (sqrt(factor3) * factor1)) * sϕ * smϕ -
            m / (factor3) * S * (ξ / factor1 * R - ∂R) * cϕ * (-cmϕ)

    #------
    mϕ = 4 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * (
        -(1 / sqrt(factor1) * ((-η / sqrt(factor3)) * ∂S + sqrt(factor3) * ∂²S) * R +
        1 / sqrt(factor3) * S * (
        sqrt(factor1) * ∂R + ξ / sqrt(factor1) * ∂²R)) * cϕ * smϕ + 
        m * (1 / sqrt(factor1) * (η / sqrt(factor3) * S +
        -1 / ((factor3)^(3 / 2)) * ∂S) * R))

    return (mξ, mη, mϕ)
end

# Nᶻ functions

function Nᶻ_oblate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    mη = 4*sqrt(factor3) / (k * d^2 * sqrt(factor2)) * (
        η * ∂S * ((2*m^2-2) * ξ / factor2^2 * R + factor1/factor2 * ∂R) -
        S * ( ((3*η^2*ξ^2 + η^2 + ξ^4 - ξ^2)/factor2^2) * ∂R + ξ * factor1 / factor2 * ∂²R) +
        m^2*ξ / (factor3 * factor1) * S * R) * cmϕ

    mξ = 4*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
    ξ * (factor3 / factor2 * ∂S + (-(2*η*(ξ^2 + 1))/(η^2 + ξ^2)^2) * S) * ∂R -
    (η * factor3 / factor2 * ∂²S + ((-η^4 - η^2 * (3*ξ^2 + 1) + ξ^2)/(η^2 + ξ^2)^2) * ∂S) * R +
    m^2*η / (factor3 * factor1) * S * R) * cmϕ

    mϕ = 4*sqrt(factor3)*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
        ξ / factor1 * ∂S * R + η / factor3 * S * ∂R) * (-smϕ)
  
    return (mξ, mη, mϕ)
end

function Nᶻ_oblate_odd(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    mη = 4*sqrt(factor3) / (k * d^2 * sqrt(factor2)) * (
        η * ∂S * ((2*m^2-2) * ξ / factor2^2 * R + factor1/factor2 * ∂R) -
        S * ( ((3*η^2*ξ^2 + η^2 + ξ^4 - ξ^2)/factor2^2) * ∂R + ξ * factor1 / factor2 * ∂²R) +
        m^2*ξ / (factor3 * factor1) * S * R) * smϕ

    mξ = 4*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
    ξ * (factor3 / factor2 * ∂S + (-(2*η*(ξ^2 + 1))/(η^2 + ξ^2)^2) * S) * ∂R -
    (η * factor3 / factor2 * ∂²S + ((-η^4 - η^2 * (3*ξ^2 + 1) + ξ^2)/(η^2 + ξ^2)^2) * ∂S) * R +
    m^2*η / (factor3 * factor1) * S * R) * smϕ

    mϕ = 4*sqrt(factor3)*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
        ξ / factor1 * ∂S * R + η / factor3 * S * ∂R) * cmϕ
  
    return (mξ, mη, mϕ)
end

function Nᶻ_prolate_even(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)
    mη = 4*sqrt(factor3) / (k * d^2 * sqrt(factor2)) * (
        η * ∂S * (-(2*m^2-2)*ξ / factor2^2 * R + factor1/factor2 * ∂R) -
        S * ( ((η^2 * (1 - 3*ξ^2) + ξ^4 + ξ^2)/(ξ^2 - η^2)^2) * ∂R + ξ * factor1 / factor2 * ∂²R) +
        m^2*ξ / (factor3 * factor1) * S * R) * cmϕ

    mξ = 4*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
    ξ * (factor3 / factor2 * ∂S + (-(2*η*(y^2 - 1))/(ξ^2 - η^2)^2) * S) * ∂R -
    (η * factor3 / factor2 * ∂²S + ((-η^4 - η^2 *(3*ξ^2 + 1) + ξ^2)/(η^2 + ξ^2)^2) * ∂S) * R +
    m^2*η / (factor3 * factor1) * S * R) * cmϕ

    mϕ = 4*sqrt(factor3)*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
        ξ / factor1 * ∂S * R + η / factor3 * S * ∂R) * (-smϕ)
  
    return (mξ, mη, mϕ)
end

function Nᶻ_prolate_odd(m, n, c, ξ, η, ϕ, k)

    d = c / k
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    sϕ, cϕ = sincos(ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    mη = 4*sqrt(factor3) / (k * d^2 * sqrt(factor2)) * (
        η * ∂S * (-(2*m^2-2)*ξ / factor2^2 * R + factor1/factor2 * ∂R) -
        S * ( ((η^2 * (1 - 3*ξ^2) + ξ^4 + ξ^2)/(ξ^2 - η^2)^2) * ∂R + ξ * factor1 / factor2 * ∂²R) +
        m^2*ξ / (factor3 * factor1) * S * R) * smϕ

    mξ = 4*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
    ξ * (factor3 / factor2 * ∂S + (-(2*η*(y^2 - 1))/(ξ^2 - η^2)^2) * S) * ∂R -
    (η * factor3 / factor2 * ∂²S + ((-η^4 - η^2 *(3*ξ^2 + 1) + ξ^2)/(η^2 + ξ^2)^2) * ∂S) * R +
    m^2*η / (factor3 * factor1) * S * R) * smϕ

    mϕ = 4*sqrt(factor3)*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
        ξ / factor1 * ∂S * R + η / factor3 * S * ∂R) * cmϕ
  
    return (mξ, mη, mϕ)
end


# Nʳ functions

function Nʳ_oblate_even(m, n, c, ξ, η, ϕ, k)

    d = c/k   
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)


    mη = 2 * sqrt(factor3) / (k * d * sqrt(factor2)) * (
        ∂S * ( ξ * factor1 / factor2 * ∂R + ((3*η^2 * ξ^2 + η^2 + ξ^4 - ξ^2)/factor2^2) * R) + 
        η * S * (factor1/factor2 * ∂²R + (2*(η^2 - 1)*ξ / factor2^2) * ∂R) - 
        m^2*η / (factor3 * factor1) * S * R) * cmϕ

        
    mξ = -2 * sqrt(factor1) / (k * d * sqrt(factor2)) * (
        +((η * factor3/factor2) *∂S + ((-η^4-η^2*(3*ξ^2+1) + ξ^2)/factor^2) * S)*∂R +
        ξ * ((-2*η * (ξ^2 + 1)/factor2) * ∂²S + factor3 / factor2 * ∂S) * R -
        m^2*ξ / (factor3 * factor1) * S * R) * cmϕ

    mϕ = 2 * m * sqrt(factor3) * sqrt(factor1) / (k * d * factor2) * (
        +1/factor1 * (η * ∂S + S) * R - 
        1 / factor3 * S * (ξ * ∂R + R)) * smϕ
    return (mξ, mη, mϕ)
end

function Nʳ_oblate_odd(m, n, c, ξ, η, ϕ, k)

    d = c/k   
    factor1 = abs(ξ)^2 + 1
    factor2 = abs(ξ)^2 + η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    mη = 2 * sqrt(factor3) / (k * d * sqrt(factor2)) * (
        ∂S * ( ξ * factor1 / factor2 * ∂R + ((3*η^2 * ξ^2 + η^2 + ξ^4 - ξ^2)/factor2^2) * R) + 
        η * S * ((factor1/factor2 * ∂²R + (2*(η^2 - 1)*ξ / factor2^2) * ∂R)) - 
        m^2*η / (factor3 * factor1) * S * R) * smϕ

    mξ = -2 * sqrt(factor1) / (k * d * sqrt(factor2)) * (
        -(((η * factor3/factor2) *∂S + ((-η^4-η^2*(3*ξ^2+1) + ξ^2)/factor^2) * S))*∂R +
        ξ * ((-2*η * (ξ^2 + 1)/factor2) * ∂²S + factor3 / factor2 * ∂S)  * R -
        m^2*ξ / (factor3 * factor1) * S * R) * smϕ

    mϕ = 2 * m * sqrt(factor3) * sqrt(factor1) / (k * d * factor2) * (
        +1/factor1 * (η * ∂S + S) * R - 
        1 / factor3 * S * (ξ * ∂R + R)) * (-cmϕ)
    return (mξ, mη, mϕ)
end

function Nʳ_prolate_even(m, n, c, ξ, η, ϕ, k)

    d = c/k   
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    mη = 2 * sqrt(factor3) / (k * d * sqrt(factor2)) * (
        ∂S * (ξ * factor1 / factor2 * ∂R + ((η^2*(1 - 3*ξ^2) + ξ^4 + ξ^2)/factor2^2) * R) - 
        η * S * (factor1/factor2 * ∂²R + (-2*factor1 * ξ / factor2^2)*∂R) + 
        m^2*η / (factor3 * factor1) * S * R) * cmϕ

    mξ = -2 * sqrt(factor1) / (k * d * sqrt(factor2)) * (
        +(((η * factor3/factor2) *∂S + ((η^4+η^2*(1 - 3*ξ^2) + ξ^2)/factor^2) * S))*∂R +
        ξ * ((-2*η * (ξ^2 - 1)/factor2) * ∂²S + factor3 / factor2 * ∂S)  * R -
        m^2*ξ / (factor3 * factor1) * S * R) * cmϕ

    mϕ = 2 * m * sqrt(factor3) * sqrt(factor1) / (k * d * factor2) * (
        -1/factor1 * (η * ∂S + S) * R - 
        1 / factor3 * S * (ξ * ∂R + R)) * smϕ
    return (mξ, mη, mϕ)
end

function Nʳ_prolate_odd(m, n, c, ξ, η, ϕ, k)

    d = c/k   
    factor1 = abs(ξ)^2 - 1
    factor2 = abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = radial_part(m, n, c, ξ, η, ϕ)
    S, ∂S, ∂²S = angular_part(m, n, c, ξ, η, ϕ)
    smϕ, cmϕ = sincos(m * ϕ)

    ξ = abs(ξ)

    mη = 2 * sqrt(factor3) / (k * d * sqrt(factor2)) * (
        ∂S * (ξ * factor1 / factor2 * ∂R + ((η^2*(1 - 3*ξ^2) + ξ^4 + ξ^2)/factor2^2) * R) - 
        η * S * (factor1/factor2 * ∂²R + (-2*factor1 * ξ / factor2^2)*∂R) + 
        m^2*η / (factor3 * factor1) * S * R) * smϕ

    mξ = -2 * sqrt(factor1) / (k * d * sqrt(factor2)) * (
        -(((η * factor3/factor2) *∂S + ((η^4+η^2*(1 - 3*ξ^2) + ξ^2)/factor^2) * S))*∂R +
        ξ * ((-2*η * (ξ^2 - 1)/factor2) * ∂²S + factor3 / factor2 * ∂S)  * R -
        m^2*ξ / (factor3 * factor1) * S * R) * smϕ

    mϕ = 2 * m * sqrt(factor3) * sqrt(factor1) / (k * d * factor2) * (
        -1/factor1 * (η * ∂S + S) * R - 
        1 / factor3 * S * (ξ * ∂R + R)) * (-cmϕ)
    return (mξ, mη, mϕ)
end