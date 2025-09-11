using TestItems
using TestItemRunner

@run_package_tests

@testsnippet Setup begin
    include("test_util.jl")
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
    EXAMPLES_DIR_ADVECTION = joinpath(examples_dir(), "linear_advection")
    EXAMPLES_DIR_BURGERS = joinpath(examples_dir(), "burgers")
    EXAMPLES_DIR_MAXWELL = joinpath(examples_dir(), "maxwell")
    EXAMPLES_DIR_EULER = joinpath(examples_dir(), "compressible_euler")
end
