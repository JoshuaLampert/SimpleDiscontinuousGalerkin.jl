using TestItems
using TestItemRunner

@run_package_tests

@testsnippet Setup begin
    include("test_util.jl")
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
end
