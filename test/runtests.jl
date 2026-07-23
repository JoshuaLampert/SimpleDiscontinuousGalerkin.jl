using TestItems
using TestItemRunner

# `SummationByPartsOperatorsExtra` requires Julia v1.11, but we also test on the Julia LTS
# (v1.10). It is therefore not a fixed dependency in `test/Project.toml` (which would make the
# test environment unresolvable on v1.10). Instead, we add it here only on v1.11+ and skip the
# tests tagged `:sbp_operators_extra` on older versions. The test environment is resolved by
# `Pkg.test` before this file runs, so this does not affect resolution on v1.10.
const SBP_OPERATORS_EXTRA_AVAILABLE = VERSION >= v"1.11"
if SBP_OPERATORS_EXTRA_AVAILABLE
    using Pkg: Pkg
    Pkg.add("SummationByPartsOperatorsExtra")
end

@run_package_tests filter=ti->SBP_OPERATORS_EXTRA_AVAILABLE ||
                              !(:sbp_operators_extra in ti.tags)

@testsnippet Setup begin
    include("test_util.jl")
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
    EXAMPLES_DIR_ADVECTION = joinpath(examples_dir(), "linear_advection")
    EXAMPLES_DIR_BURGERS = joinpath(examples_dir(), "burgers")
    EXAMPLES_DIR_MAXWELL = joinpath(examples_dir(), "maxwell")
    EXAMPLES_DIR_EULER = joinpath(examples_dir(), "compressible_euler")
end
