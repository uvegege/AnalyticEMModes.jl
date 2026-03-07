using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "..", "src", "spheroidal_utils.jl"))

@inline function _fdvec(f, x; h=1e-7, xmin=-Inf, xmax=Inf)
    xp = min(x + h, xmax)
    xm = max(x - h, xmin)
    if xp == xm
        return (f(xp) - f(x)) / max(h, eps(Float64))
    end
    return (f(xp) - f(xm)) / (xp - xm)
end

function run_spheroidal_geometry_checks(; a=1.3, Nξ=40, Nη=40, ϕ=0.63)
    zhat = SVector(0.0, 0.0, 1.0)

    function cart(type, ξ, η, ϕv)
        if type == :prolate
            x, y, z = pro2cart(a, ξ, η, ϕv)
        else
            x, y, z = obl2cart(a, ξ, η, ϕv)
        end
        return SVector(x, y, z)
    end

    function parts(type, ξ, η, ϕv)
        if type == :prolate
            vals = prolate_partials(a, ξ, η, ϕv)
            hξ, hη, hϕ = scale_factors_prolate(a, ξ, η)
            f1 = ξ^2 - 1
            f2 = ξ^2 - η^2
            f3 = 1 - η^2
        else
            vals = oblate_partials(a, ξ, η, ϕv)
            hξ, hη, hϕ = scale_factors_oblate(a, ξ, η)
            f1 = ξ^2 + 1
            f2 = ξ^2 + η^2
            f3 = 1 - η^2
        end
        rξ = SVector(vals[1], vals[2], vals[3])
        rη = SVector(vals[4], vals[5], vals[6])
        rϕ = SVector(vals[7], vals[8], vals[9])
        eξ = rξ / hξ
        eη = rη / hη
        eϕ = rϕ / hϕ
        return rξ, rη, rϕ, eξ, eη, eϕ, hξ, hη, hϕ, f1, f2, f3
    end

    function run_type(type)
        ξs = type == :prolate ? range(1.05, 3.5, length=Nξ) : range(0.05, 3.0, length=Nξ)
        ηs = range(-0.95, 0.95, length=Nη)

        max_h_from_partials = 0.0
        max_h_from_fd = 0.0
        max_orthogonality = 0.0
        max_handedness = 0.0
        max_zhat_proj = 0.0
        max_zhat_recon = 0.0
        max_ephi_z = 0.0

        for ξ in ξs, η in ηs
            rξ, rη, rϕ, eξ, eη, eϕ, hξ, hη, hϕ, f1, f2, f3 = parts(type, ξ, η, ϕ)

            # 1) scale factors vs partial norms
            max_h_from_partials = max(max_h_from_partials, abs(norm(rξ) - hξ), abs(norm(rη) - hη), abs(norm(rϕ) - hϕ))

            # 2) scale factors vs finite differences on mapping
            h = 1e-7
            ξmin = type == :prolate ? 1.0 : 0.0
            rξ_fd = _fdvec(xv -> cart(type, xv, η, ϕ), ξ; h=h, xmin=ξmin)
            rη_fd = _fdvec(yv -> cart(type, ξ, yv, ϕ), η; h=h, xmin=-1.0, xmax=1.0)
            rϕ_fd = _fdvec(pv -> cart(type, ξ, η, pv), ϕ; h=h)
            max_h_from_fd = max(max_h_from_fd, abs(norm(rξ_fd) - hξ), abs(norm(rη_fd) - hη), abs(norm(rϕ_fd) - hϕ))

            # 3) orthogonality / orientation
            max_orthogonality = max(max_orthogonality, abs(dot(eξ, eη)), abs(dot(eξ, eϕ)), abs(dot(eη, eϕ)))
            max_handedness = max(max_handedness, abs(dot(cross(eξ, eη), eϕ) - 1.0))

            # 4) z-hat components in local basis
            zξ_dot = dot(zhat, eξ)
            zη_dot = dot(zhat, eη)
            zϕ_dot = dot(zhat, eϕ)
            zξ_formula = η * sqrt(f1 / f2)
            zη_formula = ξ * sqrt(f3 / f2)
            max_zhat_proj = max(max_zhat_proj, abs(zξ_dot - zξ_formula), abs(zη_dot - zη_formula))
            max_ephi_z = max(max_ephi_z, abs(zϕ_dot))

            zrec = zξ_dot * eξ + zη_dot * eη + zϕ_dot * eϕ
            max_zhat_recon = max(max_zhat_recon, norm(zrec - zhat))
        end

        return (
            max_h_from_partials=max_h_from_partials,
            max_h_from_fd=max_h_from_fd,
            max_orthogonality=max_orthogonality,
            max_handedness=max_handedness,
            max_zhat_proj=max_zhat_proj,
            max_ephi_z=max_ephi_z,
            max_zhat_recon=max_zhat_recon,
        )
    end

    res_pro = run_type(:prolate)
    res_obl = run_type(:oblate)

    println("Geometry checks (max absolute errors)")
    println("prolate: ", res_pro)
    println("oblate : ", res_obl)

    return (prolate=res_pro, oblate=res_obl)
end

# run_spheroidal_geometry_checks()
