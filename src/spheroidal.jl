include("./spheroidal_utils.jl")

function Mᶻ(m, c, ξ, η, ϕ, basis, even, oblate)

    d = c / k
    factor1 = oblate ? abs(ξ)^2 + 1   : abs(ξ)^2 - 1
    factor2 = oblate ? abs(ξ)^2 + η^2 : abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = evaluate_radial4(basis, ξ)
    S, ∂S, ∂²S = evaluate_angular(basis, η)

    smϕ, cmϕ = sincos(m * ϕ)
    if even == false
        smϕ, cmϕ = -cmϕ, smϕ
    end

    mη = 2 * m * η / (d * sqrt(factor2) * sqrt(factor3)) * S * R * smϕ
    mξ = 2 * m * ξ / (d * sqrt(factor2) * sqrt(factor1)) * S * R * (-smϕ)
    mϕ = 2 * sqrt(factor1) * sqrt(factor3) / (d * sqrt(factor2)) * (η * ∂S * R - abs(ξ) * S * ∂R) * cmϕ

    return (mξ, mη, mϕ)
end


function Nᶻ(m, c, ξ, η, ϕ, k, basis, even, oblate)

    d = c / k
    factor1 = oblate ? abs(ξ)^2 + 1   : abs(ξ)^2 - 1
    factor2 = oblate ? abs(ξ)^2 + η^2 : abs(ξ)^2 - η^2
    factor3 = 1 - η^2

    R, ∂R, ∂²R = evaluate_radial4(basis, ξ)
    S, ∂S, ∂²S = evaluate_angular(basis, η)
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

    if even == false
        smϕ, cmϕ = -cmϕ, smϕ
    end

    ξ = abs(ξ)

    mη = 4*sqrt(factor3) / (k * d^2 * sqrt(factor2)) * (η * ∂S * (∂v1 * R + v1 * ∂R) -
        S * (∂v2 * ∂R + v2 * ∂²R) + m^2*ξ / (factor3 * factor1) * S * R) * cmϕ

    mξ = 4*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (ξ * (v3 * ∂S + ∂v3 * S) * ∂R -
    (v4 * ∂²S + ∂v4 * ∂S) * R + m^2*η / (factor3 * factor1) * S * R) * cmϕ

    mϕ = 4*sqrt(factor3)*sqrt(factor1) / (k * d^2 * sqrt(factor2)) * (
    ξ / factor1 * ∂S * R + η / factor3 * S * ∂R) * (-smϕ)
  
    return (mξ, mη, mϕ)
end
