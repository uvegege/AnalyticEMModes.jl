using Test

include(joinpath(@__DIR__, "..", "src", "spheroidal_utils.jl"))
include(joinpath(@__DIR__, "..", "src", "spheroidal_vectors.jl"))
include(joinpath(@__DIR__, "helpers", "spheroidal_comp.jl"))
include(joinpath(@__DIR__, "helpers", "spheroidal_comp2.jl"))

const _FAMILIES = (:x, :y, :z, :r)
const _CASES = (
    (m=1, n=1, c=2.4, eta=0.3, phi=0.7, even=true,  oblate=false),
    (m=1, n=1, c=4.4, eta=0.3, phi=0.7, even=false, oblate=false),
    (m=1, n=1, c=2.4, eta=0.3, phi=0.7, even=true,  oblate=true),
    (m=1, n=1, c=4.4, eta=0.3, phi=0.7, even=false, oblate=true),
)

_basis_for_case(tc) = tc.oblate ? SpheroidalB(tc.m, tc.n, complex(0.0, tc.c)) : SpheroidalB(tc.m, tc.n, tc.c)
_xis_for_case(tc) = tc.oblate ? range(0.2, 2.5, length=60) : range(1.1, 3.0, length=60)

function _assert_component_close(compres; corr_min=0.999, err_max=1e-8)
    @test real(compres.corr) >= corr_min
    @test compres.err <= err_max
end

@testset "Hardcoded vs Analitico (M families)" begin
    for tc in _CASES
        basis = _basis_for_case(tc)
        xis = _xis_for_case(tc)
        for fam in _FAMILIES
            cmp = compare_M_alpha_corr_sweep(
                fam, tc.m, tc.c, tc.eta, tc.phi, 1.0, basis;
                even=tc.even, oblate=tc.oblate, ξs=xis
            )
            _assert_component_close(cmp.ξ)
            _assert_component_close(cmp.η)
            _assert_component_close(cmp.ϕ)
        end
    end
end

@testset "Hardcoded vs Analitico (N families)" begin
    for tc in _CASES
        basis = _basis_for_case(tc)
        xis = _xis_for_case(tc)
        for fam in _FAMILIES
            cmp = compare_N_alpha_corr_sweep(
                fam, tc.m, tc.c, tc.eta, tc.phi, 1.0, basis;
                even=tc.even, oblate=tc.oblate, ξs=xis
            )
            _assert_component_close(cmp.ξ)
            _assert_component_close(cmp.η)
            _assert_component_close(cmp.ϕ)
        end
    end
end
