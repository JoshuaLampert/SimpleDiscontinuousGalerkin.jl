@testsnippet Examples begin
    examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")
end

@testitem "linear_advection.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "linear_advection.jl"),
                        l2=[0.0001574641423858981], linf=[0.00042895011704513486],
                        cons_error=[4.3021142204224816e-16])

    surface_flux = (flux_central, flux_godunov)
    @testset "Different surface flux on boundary and interior" begin
        @test_trixi_include(joinpath(examples_dir(), "linear_advection.jl"),
                            surface_flux=surface_flux,
                            l2=[0.0004272537469843011], linf=[0.0008162855376093736])
    end
end

@testitem "linear_advection_Dirichlet_boundary_condition.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_Dirichlet_boundary_condition.jl"),
                        l2=[0.00019584095949362165], linf=[0.0006856564226139367])
end

@testitem "linear_advection_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_strong_form.jl"),
                        l2=[0.00015746414238620223], linf=[0.0004289501170468002])
end

@testitem "linear_advection_FDSBP.jl" setup=[Setup, Examples] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_strong_form.jl"),
                        l2=[0.00015746414238620223], linf=[0.0004289501170468002])
end

@testitem "linear_advection_FDSBP_SAT.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl" with `surface_flux = (flux_central, flux_godunov)`
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_FDSBP_SAT.jl"),
                        boundary_conditions=boundary_condition_periodic,
                        l2=[0.00042725374698477996], linf=[0.0008162855376108169])
end

@testitem "linear_advection_FD.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_FD.jl"),
                        l2=[0.03199259748601379], linf=[0.031230470786942688])

    # Use `legendre_derivative_operator` to test if `FDSBP` gives the same result`
    coordinates_min = -1.0
    coordinates_max = 1.0
    mesh = Mesh(coordinates_min, coordinates_max, 10)
    D = legendre_derivative_operator(-1.0, 1.0, 4)
    solver = FDSBP(D, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                   volume_integral = VolumeIntegralStrongForm())
    semi = Semidiscretization(mesh, equations, initial_condition, solver)
    @testset "FDSBP with Legendre operator" begin
        # Same errors as in "linear_advection.jl"
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_FD.jl"),
                            semi=semi,
                            l2=[0.00015746414238620223], linf=[0.0004289501170468002])
    end
end

@testitem "linear_advection_flux_differencing.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_flux_differencing.jl"),
                        l2=[0.00015746414238611707], linf=[0.0004289501170472443])
end

@testitem "linear_advection_flux_differencing_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_flux_differencing_strong_form.jl"),
                        l2=[0.00015746414238603103], linf=[0.00042895011704957575])
end

@testitem "linear_advection_inhomogeneous_mesh.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_inhomogeneous_mesh.jl"),
                        l2=[0.00015505641092054123], linf=[0.0004295567744543316])

    coordinates = collect(coordinates_min:dx:coordinates_max)
    mesh = InhomogeneousMesh(coordinates)
    @testset "Inhomogeneous mesh with homogeneous nodes" begin
        @test_trixi_include(joinpath(examples_dir(),
                                     "linear_advection_inhomogeneous_mesh.jl"),
                            mesh=mesh,
                            l2=[0.0001574641423857773], linf=[0.00042895011704657815])
    end
end

@testitem "linear_advection_per_element.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                        l2=[0.0033676357851732484], linf=[0.009084487338805292],
                        cons_error=[7.494005416219807e-16])

    Ds = [isodd(element) ? D_polydeg_2 : D_polydeg_4 for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                             volume_integral = VolumeIntegralStrongForm())
    @testset "Strong form" begin
        # Same errors as above
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                            solver=solver,
                            l2=[0.0033676357851732484], linf=[0.009084487338805292])
    end

    Ds = [isodd(element) ? D_polydeg_2 : D_polydeg_4 for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralWeakForm(flux_godunov),
                             volume_integral = VolumeIntegralFluxDifferencing())
    @testset "Flux differencing" begin
        # Same errors as above
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                            solver=solver,
                            l2=[0.0033676357851732484], linf=[0.009084487338805292])
    end

    Ds = [isodd(element) ? D_polydeg_2 : D_polydeg_4 for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                             volume_integral = VolumeIntegralFluxDifferencingStrongForm())
    @testset "Flux differencing, strong form" begin
        # Same errors as above
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                            solver=solver,
                            l2=[0.0033676357851732484], linf=[0.009084487338805292])
    end

    Ds = [legendre_derivative_operator(-1.0, 1.0, 4) for element in eachelement(mesh)]
    @testset "Same operators on each element" begin
        # Same errors as in "linear_advection.jl"
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                            Ds=Ds,
                            l2=[0.0001574641423857773], linf=[0.00042895011704657815])
    end

    Ds = [legendre_derivative_operator(-1.0, 1.0, 4) for element in eachelement(mesh)]
    dx = element_spacing(mesh, 1)
    mesh = InhomogeneousMesh(collect(coordinates_min:dx:coordinates_max))
    @testset "Same operators on each element and InhomogeneousMesh" begin
        # Same errors as in "linear_advection.jl"
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_per_element.jl"),
                            Ds=Ds, mesh=mesh,
                            l2=[0.0001574641423857773], linf=[0.00042895011704657815])
    end
end
