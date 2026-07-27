function Mˣ(xi, eta, phi, S, dS, R, dR, m, d, even, type)
    if type == :oblate
        factor1 = (xi^2 + 1)
        factor2 = (xi^2 + eta^2)
    else
        factor1 = (xi^2 - 1)
        factor2 = (xi^2 - eta^2)
    end
    factor3 = (1 - eta^2)

    smϕ, cmϕ = sincos(m * phi)
    mϕ = cos(m*phi)
    if !even
        mϕ = sin(m*phi)
        smϕ, cmϕ = -cmϕ, smϕ
    end

    mx_eta = -2 * sqrt(factor1) / d / sqrt(factor2) * (S * dR * sin(phi) * cmϕ - m * xi / factor1 * S * R * cos(phi) * smϕ)
    mx_xi  =  2 * sqrt(factor3) / d / sqrt(factor2) * (dS * R * sin(phi) * cmϕ + m * eta / factor3 * S * R * cos(phi) * smϕ)
    mx_phi =  2 / (d * (factor2)) * (xi * factor3 * dS * R + eta * factor1 * S * dR) * cos(phi) * mϕ

    return mx_xi / 2, mx_eta / 2, mx_phi / 2
end


function Mʸ(xi, eta, phi, S, dS, R, dR, m, d, even, type)
    if type == :oblate
        factor1 = (xi^2 + 1)
        factor2 = (xi^2 + eta^2)
    else
        factor1 = (xi^2 - 1)
        factor2 = (xi^2 - eta^2)
    end
    factor3 = (1 - eta^2)

    smϕ, cmϕ = sincos(m * phi)
    mϕ = cos(m*phi)
    if !even
        mϕ = sin(m*phi)
        smϕ, cmϕ = -cmϕ, smϕ
    end

    my_eta = -2 * sqrt(factor1) / d / sqrt(factor2) * (-S * dR * cos(phi) * cmϕ - m * xi / factor1 * S * R * sin(phi) * smϕ)
    my_xi  =  2 * sqrt(factor3) / d / sqrt(factor2) * (-dS * R * cos(phi) * cmϕ + m * eta / factor3 * S * R * sin(phi) * smϕ)
    my_phi = 2 / (d * factor2) * (xi * factor3 * dS * R + eta * factor1 * S * dR) * sin(phi) * mϕ

    return my_xi / 2, my_eta / 2, my_phi / 2
end


function M⁺(xi, eta, phi, S, dS, R, dR, m, d, even, type)
    factor1 = type == :oblate ? xi^2 + 1 : xi^2 - 1
    factor2 = type == :oblate ? xi^2 + eta^2 : xi^2 - eta^2
    factor3 = 1 - eta^2
    sϕ, cϕ = sincos((m + 1) * phi)

    if even
        Mη = -sqrt(factor1) / (d * sqrt(factor2)) * (S * dR - m * xi / factor1 * S * R) * sϕ
        Mξ = sqrt(factor3) / (d * sqrt(factor2)) * (dS * R + m * eta / factor3 * S * R) * sϕ
        Mϕ = (xi * factor3 * dS * R + eta * factor1 * S * dR) / (d * factor2) * cϕ
    else
        Mη = sqrt(factor1) / (d * sqrt(factor2)) * (S * dR - m * xi / factor1 * S * R) * cϕ
        Mξ = -sqrt(factor3) / (d * sqrt(factor2)) * (dS * R + m * eta / factor3 * S * R) * cϕ
        Mϕ = (xi * factor3 * dS * R + eta * factor1 * S * dR) / (d * factor2) * sϕ
    end

    return Mξ / 2, Mη / 2, Mϕ / 2
end


function M⁻(xi, eta, phi, S, dS, R, dR, m, d, even, type)
    factor1 = type == :oblate ? xi^2 + 1 : xi^2 - 1
    factor2 = type == :oblate ? xi^2 + eta^2 : xi^2 - eta^2
    factor3 = 1 - eta^2
    sϕ, cϕ = sincos((m - 1) * phi)

    if even
        Mη = sqrt(factor1) / (d * sqrt(factor2)) * (S * dR + m * xi / factor1 * S * R) * sϕ
        Mξ = -sqrt(factor3) / (d * sqrt(factor2)) * (dS * R - m * eta / factor3 * S * R) * sϕ
        Mϕ = (xi * factor3 * dS * R + eta * factor1 * S * dR) / (d * factor2) * cϕ
    else
        Mη = -sqrt(factor1) / (d * sqrt(factor2)) * (S * dR + m * xi / factor1 * S * R) * cϕ
        Mξ = sqrt(factor3) / (d * sqrt(factor2)) * (dS * R - m * eta / factor3 * S * R) * cϕ
        Mϕ = (xi * factor3 * dS * R + eta * factor1 * S * dR) / (d * factor2) * sϕ
    end

    return Mξ / 2, Mη / 2, Mϕ / 2
end


function Mᶻ(xi, eta, phi, S, dS, R, dR, m, d, even, type)
    if type == :oblate
        factor1 = (xi^2 + 1)
        factor2 = (xi^2 + eta^2)
    else
        factor1 = (xi^2 - 1)
        factor2 = (xi^2 - eta^2)
    end
    factor3 = (1 - eta^2)

    smϕ, cmϕ = sincos(m * phi)
    mϕ = cos(m*phi)
    if !even
        mϕ = sin(m*phi)
        smϕ, cmϕ = -cmϕ, smϕ
    end

    mz_eta = 2 * m * eta / (d * sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mz_xi  = 2 * m * xi  / (d * sqrt(factor2) * sqrt(factor1)) * S * R * (-smϕ)
    mz_phi = 2 * sqrt(factor1) * sqrt(factor3) / (d * factor2) * (eta * dS * R - xi * S * dR) * mϕ

    return mz_xi / 2, mz_eta / 2, mz_phi / 2
end


function Mʳ(xi, eta, phi, S, dS, R, dR, m, d, even, type)
    if type == :oblate
        factor1 = (xi^2 + 1)
        factor2 = (xi^2 + eta^2)
    else
        factor1 = (xi^2 - 1)
        factor2 = (xi^2 - eta^2)
    end
    factor3 = (1 - eta^2)
    s = type == :oblate ? 1 : -1

    smϕ, cmϕ = sincos(m * phi)
    mϕ = cos(m*phi)
    if !even
        mϕ = sin(m*phi)
        smϕ, cmϕ = -cmϕ, smϕ
    end

    mr_eta = m * xi      / (sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mr_xi  = s * m * eta / (sqrt(factor2) * sqrt(factor1)) * S * R * smϕ
    mr_phi = sqrt(factor1) * sqrt(factor3) / factor2 * (xi * dS * R + s * eta * S * dR) * mϕ

    return mr_xi, mr_eta, mr_phi
end


function Nᶻ(ξ, η, ϕ, S, ∂S, ∂²S, R, ∂R, ∂²R, m, d, k, even, oblate)

    factor1 = oblate ? abs(ξ)^2 + 1   : abs(ξ)^2 - 1
    factor2 = oblate ? abs(ξ)^2 + η^2 : abs(ξ)^2 - η^2
    factor3 = 1 - η^2
    smϕ, cmϕ = sincos(m * ϕ)

    v1 = factor1 / factor2
    v2 = ξ * v1
    v3 = factor3 / factor2
    v4 = η * v3

    if oblate == true
        ∂v1 = (-2ξ*factor1) / (factor2^2) + (2ξ) / factor2
        ∂v2 = (1 + 3(ξ^2)) / factor2 + (-2*factor1*(ξ^2)) / factor2^2
        ∂v3 = (-2η) / factor2 + (-2η*factor3) / factor2^2
        ∂v4 = (1 - 3(η^2)) / factor2 + (-2*factor3*(η^2)) / factor2^2
    else
        ∂v1 = (-2ξ*factor1) / (factor2^2) + (2ξ) / factor2
        ∂v2 = (-1 + 3(ξ^2)) / factor2 + (-2*factor1*(ξ^2)) / factor2^2
        ∂v3 = (-2η) / factor2 + (2η*factor3) / factor2^2
        ∂v4 = (1 - 3(η^2)) / factor2 + (2*factor3*(η^2)) / factor2^2
    end
    
    smϕ, cmϕ = sincos(m * ϕ)
    mϕ = -sin(m*ϕ)
    if !even
        mϕ = cos(m*ϕ)
        smϕ, cmϕ = cmϕ, smϕ
    end

    ξ = abs(ξ)

    mη = 4*sqrt(factor3) / (k * d^2 * sqrt(factor2)) * (η * ∂S * (∂v1 * R + v1 * ∂R) -
        S * (∂v2 * ∂R + v2 * ∂²R) + m^2*ξ / (factor3 * factor1) * S * R) * cmϕ

    mξ = 4*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂S + ∂v3 * S) * ∂R -
    (v4 * ∂²S + ∂v4 * ∂S) * R + m^2*η / (factor3 * factor1) * S * R) * cmϕ

    mϕ = 4*m*sqrt(factor3)*sqrt(factor1) / (k * d^2 * (factor2)) * (
    ξ / factor1 * ∂S * R + η / factor3 * S * ∂R) * mϕ
  
    return (mξ / 4, mη / 4, mϕ / 4)
end


function Nˣ(ξ, η, ϕ, S, ∂S, ∂²S, R, ∂R, ∂²R, m, d, k, even, oblate)
    smϕ, cmϕ = sincos(m * ϕ)
    sm2ϕ = sin(m*ϕ)
    if !even
        sm2ϕ = -cos(m*ϕ)
        smϕ, cmϕ = cmϕ, smϕ
    end

    ξ = abs(ξ)
    η = η
    ξ2 = ξ^2
    η2 = η^2

    factor1 = oblate ? (ξ2 + 1) : (ξ2 - 1)
    factor2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    factor3 = 1 - η2
    ∂ηfactor2 = oblate ? 2η : -2η

    v1 = factor1^(3/2) / factor2
    v2 = ξ * sqrt(factor1) / factor2
    v3 = factor3^(3/2) / factor2
    v4 = η * sqrt(factor3) / factor2
    v5 = sqrt(factor3)
    v6 = sqrt(factor1)
    v7 = η / sqrt(factor3)
    v8 = ξ / sqrt(factor1)

    ∂v1 = ξ * sqrt(factor1) * (3 * factor2 - 2 * factor1) / (factor2^2)
    ∂v2 = (((factor1 + ξ2) / sqrt(factor1)) * factor2 - 2 * ξ2 * sqrt(factor1)) / (factor2^2)
    ∂v3 = -3 * η * sqrt(factor3) / factor2 - factor3^(3/2) * ∂ηfactor2 / (factor2^2)
    ∂v4 = (sqrt(factor3) - η2 / sqrt(factor3)) / factor2 - η * sqrt(factor3) * ∂ηfactor2 / (factor2^2)
    ∂v5 = -η / sqrt(factor3)
    ∂v6 = ξ / sqrt(factor1)
    ∂v7 = 1 / sqrt(factor3) + η2 / factor3^(3/2)
    ∂v8 = oblate ? (1 / factor1^(3/2)) : (-1 / factor1^(3/2))

    Nˣη = 4 / (k * d^2 * sqrt(factor2)) * ((η * S * (v1 * ∂²R + ∂v1 * ∂R) 
        - 1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (v2 * ∂R + ∂v2 * R) 
        - (m^2*η*S*R)/(factor3 * sqrt(factor1))) * cos(ϕ) * cmϕ 
        + m / sqrt(factor1) * (∂S + η/factor3 * S) * R * sin(ϕ) * sm2ϕ)

    Nˣξ = -4 / (k * d^2 * sqrt(factor2)) * ((ξ * R * (v3 * ∂²S + ∂v3 * ∂S) 
        + 1 / sqrt(factor3) * ∂R * S + factor1 * ∂R * (v4 * ∂S + ∂v4 * S) 
        - (m^2*ξ*S*R)/(factor1 * sqrt(factor3))) * cos(ϕ) * cmϕ 
        + m / sqrt(factor3) * (-∂R + ξ/factor1 * R) * S * sin(ϕ) * sm2ϕ)

    # NOTE: Li Appendix A (A-7c) appears to have a typo in the last term of N^x_ϕ:
    # the derivative must be with respect to ξ (acting on (ξ/sqrt(ξ^2 ± 1)) * R), not η.
    Nˣϕ = 4 * sqrt(factor3) * sqrt(factor1)/ (k * d^2 * (factor2)) * ((1/sqrt(factor1) * R * (v5 * ∂²S + ∂v5 * ∂S) 
        + 1 / sqrt(factor3) * S * (v6 * ∂²R + ∂v6 * ∂R)) * sin(ϕ) * cmϕ + m * ((1 / sqrt(factor1) * (v7 * ∂S + ∂v7 * S) * R 
        - 1 / sqrt(factor3) * S * (∂v8 * R + v8 * ∂R))) * cos(ϕ) * sm2ϕ)

    
    return (Nˣξ / 4, Nˣη / 4, Nˣϕ / 4)
end

function Nʸ(ξ, η, ϕ, S, ∂S, ∂²S, R, ∂R, ∂²R, m, d, k, even, oblate)
    smϕ, cmϕ = sincos(m * ϕ)
    sm2ϕ = sin(m*ϕ)
    if !even
        sm2ϕ = -cos(m*ϕ)
        smϕ, cmϕ = cmϕ, smϕ
    end

    ξ = abs(ξ)
    η = η
    ξ2 = ξ^2
    η2 = η^2

    factor1 = oblate ? (ξ2 + 1) : (ξ2 - 1)
    factor2 = oblate ? (ξ2 + η2) : (ξ2 - η2)
    factor3 = 1 - η2
    ∂ηfactor2 = oblate ? 2η : -2η

    v1 = factor1^(3/2) / factor2
    v2 = ξ * sqrt(factor1) / factor2
    v3 = factor3^(3/2) / factor2
    v4 = η * sqrt(factor3) / factor2
    v5 = sqrt(factor3)
    v6 = sqrt(factor1)
    v7 = η / sqrt(factor3)
    v8 = ξ / sqrt(factor1)

    ∂v1 = ξ * sqrt(factor1) * (3 * factor2 - 2 * factor1) / (factor2^2)
    ∂v2 = (((factor1 + ξ2) / sqrt(factor1)) * factor2 - 2 * ξ2 * sqrt(factor1)) / (factor2^2)
    ∂v3 = -3 * η * sqrt(factor3) / factor2 - factor3^(3/2) * ∂ηfactor2 / (factor2^2)
    ∂v4 = (sqrt(factor3) - η2 / sqrt(factor3)) / factor2 - η * sqrt(factor3) * ∂ηfactor2 / (factor2^2)
    ∂v5 = -η / sqrt(factor3)
    ∂v6 = ξ / sqrt(factor1)
    ∂v7 = 1 / sqrt(factor3) + η2 / factor3^(3/2)
    ∂v8 = oblate ? (1 / factor1^(3/2)) : (-1 / factor1^(3/2))

    Nʸη = 4 / (k * d^2 * sqrt(factor2)) * ((η * S * (v1 * ∂²R + ∂v1 * ∂R) 
        - 1 / sqrt(factor1) * ∂S * R + factor3 * ∂S * (v2 * ∂R + ∂v2 * R) 
        - (m^2*η*S*R)/(factor3 * sqrt(factor1))) * sin(ϕ) * cmϕ 
        - m / sqrt(factor1) * (∂S + η/factor3 * S) * R * cos(ϕ) * sm2ϕ)

    Nʸξ = -4 / (k * d^2 * sqrt(factor2)) * ((ξ * R * (v3 * ∂²S + ∂v3 * ∂S) 
        + 1 / sqrt(factor3) * ∂R * S + factor1 * ∂R * (v4 * ∂S + ∂v4 * S) 
        - (m^2*ξ*S*R)/(factor1 * sqrt(factor3))) * sin(ϕ) * cmϕ 
        - m / sqrt(factor3) * (-∂R + ξ/factor1 * R) * S * cos(ϕ) * sm2ϕ)

    # NOTE: Li Appendix A (A-7c) appears to have a typo in the last term of N^y_ϕ:
    # the derivative must be with respect to ξ (acting on (ξ/sqrt(ξ^2 ± 1)) * R), not η.
    Nʸϕ = 4 * sqrt(factor3) * sqrt(factor1)/ (k * d^2 * (factor2)) * (-(1/sqrt(factor1) * R * (v5 * ∂²S + ∂v5 * ∂S) 
        + 1 / sqrt(factor3) * S * (v6 * ∂²R + ∂v6 * ∂R)) * cos(ϕ) * cmϕ + m * ((1 / sqrt(factor1) * (v7 * ∂S + ∂v7 * S) * R 
        - 1 / sqrt(factor3) * S * (∂v8 * R + v8 * ∂R))) * sin(ϕ) * sm2ϕ)

    
    return (Nʸξ / 4, Nʸη / 4, Nʸϕ / 4)
end

function N⁺(ξ, η, ϕ, S, ∂S, ∂²S, R, ∂R, ∂²R, m, d, k, even, oblate)
    ξ = abs(ξ)
    ξ2 = ξ^2
    η2 = η^2
    factor1 = oblate ? ξ2 + 1 : ξ2 - 1
    factor2 = oblate ? ξ2 + η2 : ξ2 - η2
    factor3 = 1 - η2
    ∂ηfactor2 = oblate ? 2 * η : -2 * η

    v1 = factor1^(3/2) / factor2
    v2 = ξ * sqrt(factor1) / factor2
    v3 = factor3^(3/2) / factor2
    v4 = η * sqrt(factor3) / factor2
    v5 = sqrt(factor3)
    v6 = sqrt(factor1)
    v7 = η / sqrt(factor3)
    v8 = ξ / sqrt(factor1)

    ∂v1 = ξ * sqrt(factor1) * (3 * factor2 - 2 * factor1) / factor2^2
    ∂v2 = (((factor1 + ξ2) / sqrt(factor1)) * factor2 - 2 * ξ2 * sqrt(factor1)) / factor2^2
    ∂v3 = -3 * η * sqrt(factor3) / factor2 - factor3^(3/2) * ∂ηfactor2 / factor2^2
    ∂v4 = (sqrt(factor3) - η2 / sqrt(factor3)) / factor2 - η * sqrt(factor3) * ∂ηfactor2 / factor2^2
    ∂v5 = -η / sqrt(factor3)
    ∂v6 = ξ / sqrt(factor1)
    ∂v7 = 1 / sqrt(factor3) + η2 / factor3^(3/2)
    ∂v8 = oblate ? 1 / factor1^(3/2) : -1 / factor1^(3/2)

    sϕ, cϕ = sincos((m + 1) * ϕ)
    Nη_core = η * S * (v1 * ∂²R + ∂v1 * ∂R) + factor3 * ∂S * (v2 * ∂R + ∂v2 * R) -
              (m + 1) / sqrt(factor1) * ∂S * R - m * (m + 1) * η / (factor3 * sqrt(factor1)) * S * R
    Nξ_core = ξ * R * (v3 * ∂²S + ∂v3 * ∂S) + factor1 * ∂R * (v4 * ∂S + ∂v4 * S) +
              (m + 1) / sqrt(factor3) * S * ∂R - m * (m + 1) * ξ / (factor1 * sqrt(factor3)) * S * R
    Nϕ_core = 1 / sqrt(factor1) * R * (v5 * ∂²S + ∂v5 * ∂S) +
              1 / sqrt(factor3) * S * (v6 * ∂²R + ∂v6 * ∂R) +
              m / sqrt(factor1) * (v7 * ∂S + ∂v7 * S) * R -
              m / sqrt(factor3) * S * (∂v8 * R + v8 * ∂R)

    if even
        Nη = 2 / (k * d^2 * sqrt(factor2)) * Nη_core * cϕ
        Nξ = -2 / (k * d^2 * sqrt(factor2)) * Nξ_core * cϕ
        Nϕ = 2 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * Nϕ_core * sϕ
    else
        Nη = 2 / (k * d^2 * sqrt(factor2)) * Nη_core * sϕ
        Nξ = -2 / (k * d^2 * sqrt(factor2)) * Nξ_core * sϕ
        Nϕ = -2 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * Nϕ_core * cϕ
    end

    return Nξ / 4, Nη / 4, Nϕ / 4
end

function N⁻(ξ, η, ϕ, S, ∂S, ∂²S, R, ∂R, ∂²R, m, d, k, even, oblate)
    ξ = abs(ξ)
    ξ2 = ξ^2
    η2 = η^2
    factor1 = oblate ? ξ2 + 1 : ξ2 - 1
    factor2 = oblate ? ξ2 + η2 : ξ2 - η2
    factor3 = 1 - η2
    ∂ηfactor2 = oblate ? 2 * η : -2 * η

    v1 = factor1^(3/2) / factor2
    v2 = ξ * sqrt(factor1) / factor2
    v3 = factor3^(3/2) / factor2
    v4 = η * sqrt(factor3) / factor2
    v5 = sqrt(factor3)
    v6 = sqrt(factor1)
    v7 = η / sqrt(factor3)
    v8 = ξ / sqrt(factor1)

    ∂v1 = ξ * sqrt(factor1) * (3 * factor2 - 2 * factor1) / factor2^2
    ∂v2 = (((factor1 + ξ2) / sqrt(factor1)) * factor2 - 2 * ξ2 * sqrt(factor1)) / factor2^2
    ∂v3 = -3 * η * sqrt(factor3) / factor2 - factor3^(3/2) * ∂ηfactor2 / factor2^2
    ∂v4 = (sqrt(factor3) - η2 / sqrt(factor3)) / factor2 - η * sqrt(factor3) * ∂ηfactor2 / factor2^2
    ∂v5 = -η / sqrt(factor3)
    ∂v6 = ξ / sqrt(factor1)
    ∂v7 = 1 / sqrt(factor3) + η2 / factor3^(3/2)
    ∂v8 = oblate ? 1 / factor1^(3/2) : -1 / factor1^(3/2)

    sϕ, cϕ = sincos((m - 1) * ϕ)
    Nη_core = η * S * (v1 * ∂²R + ∂v1 * ∂R) + factor3 * ∂S * (v2 * ∂R + ∂v2 * R) +
              (m - 1) / sqrt(factor1) * ∂S * R - m * (m - 1) * η / (factor3 * sqrt(factor1)) * S * R
    Nξ_core = ξ * R * (v3 * ∂²S + ∂v3 * ∂S) + factor1 * ∂R * (v4 * ∂S + ∂v4 * S) -
              (m - 1) / sqrt(factor3) * S * ∂R - m * (m - 1) * ξ / (factor1 * sqrt(factor3)) * S * R
    Nϕ_core = 1 / sqrt(factor1) * R * (v5 * ∂²S + ∂v5 * ∂S) +
              1 / sqrt(factor3) * S * (v6 * ∂²R + ∂v6 * ∂R) -
              m / sqrt(factor1) * (v7 * ∂S + ∂v7 * S) * R +
              m / sqrt(factor3) * S * (∂v8 * R + v8 * ∂R)

    if even
        Nη = 2 / (k * d^2 * sqrt(factor2)) * Nη_core * cϕ
        Nξ = -2 / (k * d^2 * sqrt(factor2)) * Nξ_core * cϕ
        Nϕ = -2 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * Nϕ_core * sϕ
    else
        Nη = 2 / (k * d^2 * sqrt(factor2)) * Nη_core * sϕ
        Nξ = -2 / (k * d^2 * sqrt(factor2)) * Nξ_core * sϕ
        Nϕ = 2 * sqrt(factor3) * sqrt(factor1) / (k * d^2 * factor2) * Nϕ_core * cϕ
    end

    return Nξ / 4, Nη / 4, Nϕ / 4
end

function Nʳ(ξ, η, ϕ, S, ∂S, ∂²S, R, ∂R, ∂²R, m, d, k, even, oblate)
    smϕ, cmϕ = sincos(m * ϕ)
    sm2ϕ = sin(m * ϕ)
    if !even
        sm2ϕ = -cos(m * ϕ)
        smϕ, cmϕ = cmϕ, smϕ
    end

    ξ = abs(ξ)
    ξ2 = ξ^2
    η2 = η^2

    factor1 = oblate ? (ξ2 + 1) : (ξ2 - 1)  
    factor2 = oblate ? (ξ2 + η2) : (ξ2 - η2) 
    factor3 = 1 - η2
    sgn = oblate ? 1.0 : -1.0  
    ∂ηfactor2 = oblate ? 2η : -2η

    v1 = ξ * factor1 / factor2
    v2 = factor1 / factor2
    ∂v1 = ((factor1 + 2ξ2) * factor2 - 2ξ2 * factor1) / factor2^2
    ∂v2 = 2ξ * (factor2 - factor1) / factor2^2
    v3 = η * factor3 / factor2
    v4 = factor3 / factor2
    ∂v3 = ((factor3 - 2η2) * factor2 - η * factor3 * ∂ηfactor2) / factor2^2
    ∂v4 = (-2η * factor2 - factor3 * ∂ηfactor2) / factor2^2

    Nʳη = 2 * sqrt(factor3) / (k * d * sqrt(factor2)) * (
        ∂S * (∂v1 * R + v1 * ∂R) +
        sgn * η * S * (∂v2 * ∂R + v2 * ∂²R) -
        sgn * (m^2 * η / (factor3 * factor1)) * S * R) * cmϕ

    Nʳξ = -2 * sqrt(factor1) / (k * d * sqrt(factor2)) * (
        sgn * (∂v3 * S + v3 * ∂S) * ∂R +
        ξ * (∂v4 * ∂S + v4 * ∂²S) * R -
        (m^2 * ξ / (factor3 * factor1)) * S * R) * cmϕ

    Nʳϕ = 2 * m * sqrt(factor3) * sqrt(factor1) / (k * d * factor2) * (
        sgn * ((S + η * ∂S) * R) / factor1 -
        (S * (R + ξ * ∂R)) / factor3) * sm2ϕ

    return (Nʳξ / 2, Nʳη / 2, Nʳϕ / 2)
end
