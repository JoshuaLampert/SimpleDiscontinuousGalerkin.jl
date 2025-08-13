@testitem "utils" begin
    # Interpolation
    D = legendre_derivative_operator(-1.0, 1.0, 5)
    x = 0.5
    nodes = grid(D)
    values = sin.(nodes)
    baryweights = PolynomialBases.barycentric_weights(nodes)
    f = interpolate(x, values, nodes, baryweights)

    e_M = SimpleDiscontinuousGalerkin.interpolation_operator(x, D)
    @test_nowarn isapprox(e_M' * values, sin(x), atol = 1.0e-14)
    @test_nowarn isapprox(e_M' * values, f, atol = 1.0e-14)

    D = derivative_operator(MattssonNordström2004(), 1, 2, -1.1, 2.6, 20)
    nodes = grid(D)
    values = sin.(nodes)

    e_M = SimpleDiscontinuousGalerkin.interpolation_operator(x, D)
    @test_nowarn isapprox(e_M' * values, sin(x), atol = 1.0e-14)

    # convergence test
    convergence_test(default_example(), 3, interval = 1000)
    polydegs = [1, 3, 5]
    for polydeg in polydegs
        eoc_mean_values, _ = convergence_test(default_example(), 3, N_elements = 16,
                                              tspan = (0.0, 1.0), polydeg = polydeg,
                                              abstol = 1e-14, reltol = 1e-14,
                                              interval = 1000)
        @test isapprox(eoc_mean_values[:l2][1], polydeg + 1, atol = 0.15)
        @test isapprox(eoc_mean_values[:linf][1], polydeg + 1, atol = 0.15)

        eoc_mean_values2, _ = convergence_test(default_example(), [16, 32, 64],
                                               tspan = (0.0, 1.0), polydeg = polydeg,
                                               abstol = 1e-14, reltol = 1e-14,
                                               interval = 1000)
        for kind in (:l2, :linf), variable in (1,)
            eoc_mean_values[kind][variable] == eoc_mean_values2[kind][variable]
        end
    end
end

@testitem "equations" begin
    equations = @test_nowarn LinearAdvectionEquation1D(-2.0)
    @test_nowarn print(equations)
    @test_nowarn display(equations)
    @test ndims(equations) == 1
    @test SimpleDiscontinuousGalerkin.nvariables(equations) == 1
    @test SimpleDiscontinuousGalerkin.get_name(equations) == "LinearAdvectionEquation1D"

    u1 = [SVector(1.0), SVector(-2.0)]
    for conversion in (cons2cons, cons2prim)
        @test length(varnames(conversion, equations)) ==
              length(conversion(first(u1), equations))
    end
    @test cons2cons.(u1, equations) == u1
    @test cons2entropy.(u1, equations) == u1
    @test prim2cons.(cons2prim.(u1, equations), equations) == u1
    @test flux.(u1, equations) == [SVector(-2.0), SVector(4.0)]
    @test flux_central.(u1, u1, equations) == flux.(u1, equations)
    @test flux_godunov.(u1, u1, equations) == flux.(u1, equations)
    @test flux_lax_friedrichs.(u1, u1, equations) == flux.(u1, equations)

    equations = @test_nowarn BurgersEquation1D()
    for conversion in (cons2cons, cons2prim)
        @test length(varnames(conversion, equations)) ==
              length(conversion(first(u1), equations))
    end
    @test cons2cons.(u1, equations) == u1
    @test cons2entropy.(u1, equations) == u1
    @test prim2cons.(cons2prim.(u1, equations), equations) == u1
    @test flux_godunov.(u1, u1, equations) == flux.(u1, equations)
    @test flux_lax_friedrichs.(u1, u1, equations) == flux.(u1, equations)
    @test flux_ec.(u1, u1, equations) == flux.(u1, equations)

    u2 = [SVector(1.0, 2.0), SVector(-2.0, -3.0)]
    equations = @test_nowarn MaxwellEquations1D(3.12)
    for conversion in (cons2cons, cons2prim)
        @test length(varnames(conversion, equations)) ==
              length(conversion(first(u2), equations))
    end
    @test cons2cons.(u2, equations) == u2
    @test all(isapprox.(cons2entropy.(u2, equations),
                        [SVector(1.0, 19.4688), SVector(-2.0, -29.2032)]))
    @test prim2cons.(cons2prim.(u2, equations), equations) == u2
    @test electric_field.(u2, equations) == [1.0, -2.0]
    @test magnetic_field.(u2, equations) == [2.0, -3.0]
    @test flux_godunov.(u2, u2, equations) == flux.(u2, equations)
    @test flux_lax_friedrichs.(u2, u2, equations) == flux.(u2, equations)

    u3 = [SVector(1.0, 2.0, 4.0), SVector(2.0, 3.0, 3.0)]
    equations = @test_nowarn CompressibleEulerEquations1D(1.4)
    for conversion in (cons2cons, cons2prim)
        @test length(varnames(conversion, equations)) ==
              length(conversion(first(u3), equations))
    end
    @test cons2cons.(u3, equations) == u3
    @test all(isapprox.(cons2entropy.(u3, equations),
                        [
                            SVector(1.5578588782855252, 2.5, -1.25),
                            SVector(1.4359471427746477, 10.0, -6.666666666666668)
                        ]))
    @test prim2cons.(cons2prim.(u3, equations), equations) == u3
    @test density.(u3, equations) == [1.0, 2.0]
    @test velocity.(u3, equations) == [2.0, 1.5]
    @test momentum.(u3, equations) == [2.0, 3.0]
    @test all(isapprox.(pressure.(u3, equations), [0.8, 0.3]))
    @test all(isapprox.(density_pressure.(u3, equations), [0.8, 0.6]))
    @test all(isapprox.(entropy_thermodynamic.(u3, equations),
                        [-0.22314355131421, -2.1743788571098595]))
    @test all(isapprox.(entropy_math.(u3, equations),
                        [0.557858878285525, 10.871894285549299]))
    @test energy_total.(u3, equations) == [4.0, 3.0]
    @test all(isapprox.(flux_lax_friedrichs.(u3, u3, equations), flux.(u3, equations)))
    @test all(isapprox.(flux_hll.(u3, u3, equations), flux.(u3, equations)))
    @test all(isapprox.(flux_ranocha.(u3, u3, equations), flux.(u3, equations)))
    @test all(isapprox.(flux_godunov.(u3, u3, equations), flux.(u3, equations)))

    @test_nowarn print(FluxLaxFriedrichs())
    @test_nowarn display(FluxLaxFriedrichs())
    @test_nowarn print(FluxHLL())
    @test_nowarn display(FluxHLL())
    @test_nowarn print(SimpleDiscontinuousGalerkin.DissipationLocalLaxFriedrichs())
    @test_nowarn display(SimpleDiscontinuousGalerkin.DissipationLocalLaxFriedrichs())
    @test_nowarn print(SimpleDiscontinuousGalerkin.FluxPlusDissipation(flux_godunov,
                                                                       SimpleDiscontinuousGalerkin.DissipationLocalLaxFriedrichs()))
    @test_nowarn display(SimpleDiscontinuousGalerkin.FluxPlusDissipation(flux_godunov,
                                                                         SimpleDiscontinuousGalerkin.DissipationLocalLaxFriedrichs()))
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
    @test_nowarn show(IOContext(stdout, :compact => false), mesh)
    @test ndims(mesh) == 1
    @test SimpleDiscontinuousGalerkin.nelements(mesh) == N
    @test real(mesh) == Float64
    @test SimpleDiscontinuousGalerkin.xmin(mesh) == coordinates_min
    @test SimpleDiscontinuousGalerkin.xmax(mesh) == coordinates_max

    mesh = InhomogeneousMesh(coordinates_min:0.1:coordinates_max)
    @test_nowarn print(mesh)
    @test_nowarn display(mesh)
    @test_nowarn show(IOContext(stdout, :compact => false), mesh)
    @test ndims(mesh) == 1
    @test SimpleDiscontinuousGalerkin.nelements(mesh) == 10
    @test real(mesh) == Float64
    @test SimpleDiscontinuousGalerkin.xmin(mesh) == coordinates_min
    @test SimpleDiscontinuousGalerkin.xmax(mesh) == coordinates_max

    mesh_left = Mesh(-1.0, 0.1, 5)
    mesh_right = Mesh(-0.1, 1.0, 5)
    mesh = OversetGridMesh(mesh_left, mesh_right)
    @test_nowarn print(mesh)
    @test_nowarn display(mesh)
    @test_nowarn show(IOContext(stdout, :compact => false), mesh)
    @test ndims(mesh) == 1
    @test SimpleDiscontinuousGalerkin.nelements(mesh) == 10
    @test real(mesh) == Float64
    @test SimpleDiscontinuousGalerkin.xmin(mesh) == -1.0
    @test SimpleDiscontinuousGalerkin.xmax(mesh) == 1.0
    @test SimpleDiscontinuousGalerkin.left_overlap_element(mesh) == 5
    @test SimpleDiscontinuousGalerkin.right_overlap_element(mesh) == 1
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
    @test real(semi) == Float64
    @test SimpleDiscontinuousGalerkin.nvariables(semi) == 1
    @test SimpleDiscontinuousGalerkin.nelements(semi) == 5
    @test SimpleDiscontinuousGalerkin.ndofs(semi) == 20
    @test nnodes(semi, 2) == 4
    @test eachnode(semi, 2) == 1:4

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

@testitem "Jacobian" begin
    using LinearAlgebra: eigvals
    trixi_include(@__MODULE__, joinpath(examples_dir(), "linear_advection.jl"),
                  tspan = (0.0, 0.01))
    J = @test_nowarn jacobian_fd(semi)
    # This is stable
    @test maximum(real, eigvals(J)) < 0.0

    trixi_include(@__MODULE__, joinpath(examples_dir(), "linear_advection.jl"),
                  tspan = (0.0, 0.01), surface_flux = flux_central)
    J = @test_nowarn jacobian_fd(semi)
    # This is conservative
    @test maximum(abs.(real.(eigvals(J)))) < 1e-7

    trixi_include(@__MODULE__, joinpath(examples_dir(), "linear_advection_per_element.jl"),
                  tspan = (0.0, 0.01))
    J = @test_nowarn jacobian_fd(semi)
    # This is stable
    @test maximum(real, eigvals(J)) < 1e-7

    trixi_include(@__MODULE__, joinpath(examples_dir(), "linear_advection_overset_grid.jl"),
                  tspan = (0.0, 0.01))
    J = @test_nowarn jacobian_fd(semi)
    # This has some eigenvalues with slightly positive real part
    @test count(real.(eigvals(J)) .> 1e-7) == 10
    @test maximum(real, eigvals(J)) < 1e-3

    trixi_include(@__MODULE__, joinpath(examples_dir(), "maxwell_overset_grid.jl"),
                  tspan = (0.0, 0.01))
    J = @test_nowarn jacobian_fd(semi)
    # This has some eigenvalues with slightly positive real part
    @test count(real.(eigvals(J)) .> 1e-7) == 20
    @test maximum(real, eigvals(J)) < 1e-3
end

@testitem "SummaryCallback" begin
    summary_callback = SummaryCallback()
    @test_nowarn print(summary_callback)
    @test_nowarn display(summary_callback)
end

@testitem "AnalysisCallback" setup=[Setup] begin
    trixi_include(@__MODULE__, joinpath(examples_dir(), "linear_advection.jl"),
                  tspan = (0.0, 0.01))
    @test_nowarn print(analysis_callback)
    @test_nowarn display(analysis_callback)
    l2, linf = analysis_callback(sol)
    errs = errors(analysis_callback)
    @test l2[1] == errs.l2_error[end]
    @test linf[1] == errs.linf_error[end]
end

@testitem "StepsizeCallback" setup=[Setup] begin
    stepsize_callback = StepsizeCallback(cfl = 1.0)
    @test_nowarn print(stepsize_callback)
    @test_nowarn display(stepsize_callback)
end

@testitem "visualization" setup=[Setup] begin
    using Plots
    include(joinpath(examples_dir(), "linear_advection.jl"))
    @test_nowarn plot(flat_grid(semi), get_variable(sol.u[end], 1, semi))
    @test_nowarn plot(semi => sol)
    @test_nowarn plot(semi => sol, plot_initial = true)
    @test_nowarn plot(semi => sol, step = 5)
    @test_nowarn plot(semi, sol, plot_initial = true, step = 6)

    include(joinpath(examples_dir(), "linear_advection_per_element.jl"))
    @test_nowarn plot(flat_grid(semi), get_variable(sol.u[end], 1, semi))
    @test_nowarn plot(semi, sol, plot_initial = true, step = 6)
    @test_nowarn plot(analysis_callback)
    @test_nowarn plot(analysis_callback, what = (:errors,))
    @test_nowarn plot(analysis_callback, what = (:integrals, :errors))

    include(joinpath(examples_dir(), "linear_advection_overset_grid.jl"))
    @test_nowarn plot(semi => sol, plot_initial = true, step = 6)

    include(joinpath(examples_dir(), "linear_advection_overset_grid_per_element.jl"))
    @test_nowarn plot(semi => sol, plot_initial = true, step = 6)
end
