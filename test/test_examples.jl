@testsnippet Examples begin
    using TrixiBase: trixi_include
    using TrixiTest: @trixi_test_nowarn
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
end

@testitem "linear_advection.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(), "linear_advection.jl"))
end
