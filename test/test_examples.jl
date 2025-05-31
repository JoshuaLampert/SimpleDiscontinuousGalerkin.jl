@testsnippet Examples begin
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
end

@testitem "linear_advection.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "linear_advection.jl"),
                        l2=[0.00016123961626763236], linf=[0.00042895011704513486])

    surface_flux = (flux_central, flux_godunov)
    @test_trixi_include(joinpath(examples_dir(), "linear_advection.jl"),
                        surface_flux=surface_flux,
                        l2=[0.00041060056908958845], linf=[0.0008162855376093736])
end

@testitem "linear_advection_Dirichlet_boundary_condition.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_Dirichlet_boundary_condition.jl"),
                        l2=[0.00015725543332057758], linf=[0.0006856564226139367])
end

@testitem "linear_advection_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_strong_form.jl"),
                        l2=[0.00016123961626768018], linf=[0.0004289501170468002])
end

@testitem "linear_advection_FDSBP.jl" setup=[Setup, Examples] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_strong_form.jl"),
                        l2=[0.00016123961626768018], linf=[0.0004289501170468002])
end

@testitem "linear_advection_FDSBP_SAT.jl" setup=[Setup] begin
    # (Almost) same errors as in "linear_advection.jl" with `surface_flux = (flux_central, flux_godunov)`
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_FDSBP_SAT.jl"),
                        l2=[0.00040873810317442036], linf=[0.0008088018130785191])
end

@testitem "linear_advection_FD.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_FD.jl"),
                        l2=[0.015581091289174379], linf=[0.031230470786942688])

    # Use `legendre_derivative_operator` to test if `FDSBP` gives the same result`
    coordinates_min = -1.0
    coordinates_max = 1.0
    mesh = Mesh(coordinates_min, coordinates_max, 10)
    D = legendre_derivative_operator(-1.0, 1.0, 4)
    solver = FDSBP(D, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                   volume_integral = VolumeIntegralStrongForm())
    semi = Semidiscretization(mesh, equations, initial_condition, solver)
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_FD.jl"),
                        semi=semi,
                        l2=[0.00016123961626768018], linf=[0.0004289501170468002])
end

@testitem "linear_advection_flux_differencing.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_flux_differencing.jl"),
                        l2=[0.00016123961626813603], linf=[0.0004289501170472443])
end

@testitem "linear_advection_flux_differencing_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_flux_differencing_strong_form.jl"),
                        l2=[0.00016123961626769392], linf=[0.00042895011704957575])
end

@testitem "linear_advection_per_element.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                        l2=[0.23406659150246234], linf=[0.564353999854611])

    Ds = [legendre_derivative_operator(-1.0, 1.0, 4) for element in eachelement(mesh)]
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                        Ds=Ds,
                        l2=[0.00016123961626746092], linf=[0.00042895011704657815])
end
