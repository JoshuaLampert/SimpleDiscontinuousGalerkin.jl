@testitem "equations" begin
    equations = @test_nowarn LinearAdvectionEquation1D(-2.0)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    @test ndims(equations) == 1
    @test SimpleDiscontinuousGalerkin.nvariables(equations) == 1
    @test SimpleDiscontinuousGalerkin.get_name(equations) == "LinearAdvectionEquation1D"
    u = [SVector(1.0), SVector(-2.0)]
    @test cons2cons.(u, equations) == u
    @test flux.(u, equations) == [SVector(2.0), SVector(-4.0)]
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
    @test_nowarn print(solver)
    @test_nowarn display(solver)

    D = derivative_operator(MattssonNordström2004(), 1, 2, -1.0, 1.0, 4)
    solver = FDSBP(D)
    @test_nowarn print(solver)
    @test_nowarn display(solver)
end

@testitem "semidiscretization" begin
    equations = LinearAdvectionEquation1D(2.0)
    mesh = Mesh(0.0, 1.0, 10)
    boundary_conditions = (boundary_condition_do_nothing, boundary_condition_periodic)
    solver = DGSEM(polydeg = 3)
    semi = Semidiscretization(mesh, equations, initial_condition_convergence_test, solver;
                              boundary_conditions)
    @test_nowarn print(semi)
    @test_nowarn display(semi)
    @test ndims(semi) == 1
    @test SimpleDiscontinuousGalerkin.nvariables(semi) == 1
    @test SimpleDiscontinuousGalerkin.nelements(semi) == 10
    @test SimpleDiscontinuousGalerkin.ndofs(semi) == 40
    @test real(semi) == Float64
end
