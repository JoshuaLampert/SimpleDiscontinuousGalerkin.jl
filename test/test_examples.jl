@testsnippet Examples begin
    using TrixiBase: trixi_include
    using TrixiTest: @trixi_test_nowarn
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
end

@testitem "linear_advection.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(), "linear_advection.jl"))
end

@testitem "linear_advection_Dirichlet_boundary_condition.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_Dirichlet_boundary_condition.jl"))
end

@testitem "linear_advection_strong_form.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_strong_form.jl"))
end

@testitem "linear_advection_FDSBP.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_FDSBP.jl"))
end

@testitem "linear_advection_FDSBP_SAT.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_FDSBP_SAT.jl"))
end

@testitem "linear_advection_FD.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_FD.jl"))
end

@testitem "linear_advection_flux_differencing.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_flux_differencing.jl"))
end

@testitem "linear_advection_flux_differencing_strong_form.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_flux_differencing_strong_form.jl"))
end

@testitem "linear_advection_per_element.jl" setup=[Examples] begin
    @trixi_test_nowarn trixi_include(joinpath(examples_dir(),
                                              "linear_advection_per_element.jl"))
end
