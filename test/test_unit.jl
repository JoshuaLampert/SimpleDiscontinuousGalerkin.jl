@testitem "equations" begin
    equations = @test_nowarn LinearAdvectionEquation1D(-2.0)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    @test ndims(equations) == 1
    @test SimpleDiscontinuousGalerkin.nvariables(equations) == 1
    @test SimpleDiscontinuousGalerkin.get_name(equations) == "LinearAdvectionEquation1D"
    u = [SVector(1.0), SVector(-2.0)]
    @test cons2cons.(u, equations) == u
    @test flux.(u, equations) == [SVector(-2.0), SVector(4.0)]
    @test flux_central.(u, u, equations) == flux.(u, equations)
    @test flux_godunov.(u, u, equations) == flux.(u, equations)
end

@testitem "boundary_conditions" begin
    @test_nowarn print(boundary_condition_periodic)
    @test_nowarn print(boundary_condition_do_nothing)
    @test_nowarn print(BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0)))

    @test SimpleDiscontinuousGalerkin.digest_boundary_conditions((boundary_condition_do_nothing,
                                                                  boundary_condition_periodic)) ==
          SimpleDiscontinuousGalerkin.digest_boundary_conditions((x_neg = boundary_condition_do_nothing,
                                                                  x_pos = boundary_condition_periodic))
end

@testitem "mesh" begin
    coordinates_min = 0.0
    coordinates_max = 1.0
    N = 10
    mesh = Mesh(coordinates_min, coordinates_max, N)
    @test_nowarn print(mesh)
    @test_nowarn display(mesh)
    @test ndims(mesh) == 1
    @test SimpleDiscontinuousGalerkin.nelements(mesh) == N
    @test real(mesh) == Float64
    @test SimpleDiscontinuousGalerkin.xmin(mesh) == coordinates_min
    @test SimpleDiscontinuousGalerkin.xmax(mesh) == coordinates_max
end

@testitem "surface integrals" begin
    @test SurfaceIntegralStrongForm((flux_central, flux_central)) ==
          SurfaceIntegralStrongForm(flux_central) == SurfaceIntegralStrongForm()
    integral = SurfaceIntegralStrongForm(flux_central)
    @test_nowarn print(integral)
    @test_nowarn display(integral)

    @test SurfaceIntegralWeakForm((flux_central, flux_central)) ==
          SurfaceIntegralWeakForm(flux_central) == SurfaceIntegralWeakForm()
    integral = SurfaceIntegralWeakForm(flux_central)
    @test_nowarn print(integral)
    @test_nowarn display(integral)
end

@testitem "volume integrals" begin
    for integral in (VolumeIntegralStrongForm(),
                     VolumeIntegralWeakForm(),
                     VolumeIntegralFluxDifferencing(flux_central),
                     VolumeIntegralFluxDifferencingStrongForm(flux_central))
        @test_nowarn print(integral)
        @test_nowarn display(integral)
    end
end

@testitem "solvers" begin
    solver = DGSEM(polydeg = 3)
    @test_nowarn summary(solver)
    @test_nowarn print(solver)
    @test_nowarn display(solver)

    D = derivative_operator(MattssonNordström2004(), 1, 2, -1.0, 1.0, 4)
    solver = FDSBP(D)
    @test_nowarn summary(solver)
    @test_nowarn print(solver)
    @test_nowarn display(solver)

    Ds = [derivative_operator(MattssonNordström2004(), 1, 2, -1.0, 1.0, 4),
        derivative_operator(MattssonNordström2004(), 1, 2, -1.0, 1.0, 5)]
    solver = PerElementFDSBP(Ds)
    @test_nowarn summary(solver)
    @test_nowarn print(solver)
    @test_nowarn display(solver)
end

@testitem "semidiscretization" begin
    equations = LinearAdvectionEquation1D(2.0)
    mesh = Mesh(0.0, 1.0, 5)
    boundary_conditions = (boundary_condition_do_nothing, boundary_condition_periodic)
    D_leg = legendre_derivative_operator(-2.0, 3.0, 4)
    solver = DGSEM(polydeg = 3)
    semi = Semidiscretization(mesh, equations, initial_condition_convergence_test, solver;
                              boundary_conditions)
    @test_nowarn print(semi)
    @test_nowarn display(semi)
    @test ndims(semi) == 1
    @test SimpleDiscontinuousGalerkin.nvariables(semi) == 1
    @test SimpleDiscontinuousGalerkin.nelements(semi) == 5
    @test SimpleDiscontinuousGalerkin.ndofs(semi) == 20
    @test real(semi) == Float64

    @test all(isapprox.(grid(semi),
                        [0.0 0.2 0.4 0.6 0.8;
                         0.05527864045000422 0.25527864045000426 0.4552786404500042 0.6552786404500043 0.8552786404500042;
                         0.14472135954999582 0.34472135954999583 0.5447213595499958 0.7447213595499959 0.9447213595499958;
                         0.2 0.4 0.6 0.8 1.0], atol = 1.0e-14))

    using RecursiveArrayTools: VectorOfArray
    D_leg_2 = legendre_derivative_operator(-1.0, 1.0, 3)
    Ds = [isodd(element) ? D_leg_2 : D_leg for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_flux = flux_godunov)
    semi = Semidiscretization(mesh, equations, initial_condition_convergence_test, solver;
                              boundary_conditions)
    @test isapprox(grid(semi),
                   VectorOfArray([[0.0, 0.1, 0.2],
                                     [0.2, 0.2552786404500042, 0.34472135954999583, 0.4],
                                     [0.4, 0.5, 0.6],
                                     [0.6, 0.6552786404500043, 0.7447213595499959, 0.8],
                                     [0.8, 0.9, 1.0]]), atol = 1.0e-14)
end

@testitem "SummaryCallback" begin
    summary_callback = SummaryCallback()
    @test_nowarn print(summary_callback)
    @test_nowarn display(summary_callback)
end

@testitem "AnalysisCallback" begin
    semi = Semidiscretization(Mesh(0.0, 1.0, 5), LinearAdvectionEquation1D(2.0),
                              initial_condition_convergence_test, DGSEM(polydeg = 3))
    analysis_callback = AnalysisCallback(semi; interval = 10,
                                         extra_analysis_integrals = (mass, entropy))
    @test_nowarn print(analysis_callback)
    @test_nowarn display(analysis_callback)
end

@testitem "visualization" setup=[Setup] begin
    using Plots
    include(joinpath(examples_dir(), "linear_advection.jl"))
    @test_nowarn plot(flat_grid(semi), get_variable(sol.u[end], 1, semi))
    include(joinpath(examples_dir(), "linear_advection_per_element.jl"))
    @test_nowarn plot(flat_grid(semi), get_variable(sol.u[end], 1, semi))
    @test_nowarn plot(analysis_callback)
    @test_nowarn plot(analysis_callback, what = (:errors,))
    @test_nowarn plot(analysis_callback, what = (:integrals, :errors))
end
