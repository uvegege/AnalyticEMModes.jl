using Test

include("basic_exports.jl")

if get(ENV, "AEM_RUN_LOCAL_TESTS", "0") == "1"
    include("spheroidal_hardcoded_vs_analytic_test.jl")
end
