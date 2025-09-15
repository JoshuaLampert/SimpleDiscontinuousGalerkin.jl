@testitem "linear_advection.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection.jl"),
                        l2=[0.0001574641423858981], linf=[0.00042895011704513486],
                        cons_error=[4.3021142204224816e-16],
                        change_mass=-3.0531133177191805e-16,
                        change_entropy=-9.082223088596741e-7,
                        entropy_timederivative=-9.210638494683288e-7)

    surface_flux = (flux_central, flux_godunov)
    @testset "Different surface flux on boundary and interior" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection.jl"),
                            surface_flux=surface_flux,
                            l2=[0.0004272537469843011], linf=[0.0008162855376093736],
                            cons_error=[7.7021722333370235e-16],
                            change_mass=7.7021722333370235e-16,
                            change_entropy=-4.29859495665319e-7,
                            entropy_timederivative=-6.420623105407586e-7)
    end

    surface_flux = flux_lax_friedrichs
    @testset "flux_lax_friedrichs as surface flux" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection.jl"),
                            surface_flux=surface_flux,
                            l2=[0.00015746414238601946], linf=[0.00042895011704702224],
                            cons_error=[3.191891195797325e-16],
                            change_mass=3.191891195797325e-16,
                            change_entropy=-9.082223089151853e-7,
                            entropy_timederivative=-9.210638497458845e-7)
    end
end

@testitem "linear_advection_Dirichlet_boundary_condition.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                 "linear_advection_Dirichlet_boundary_condition.jl"),
                        l2=[0.00019584095949362165], linf=[0.0006856564226139367],
                        change_mass=-0.5462956322467551,
                        change_entropy=-0.1978580355978945,
                        entropy_timederivative=-0.0066261099903804676)
end

@testitem "linear_advection_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_strong_form.jl"),
                        l2=[0.00015746414238620223], linf=[0.0004289501170468002],
                        change_mass=2.7755575615628914e-17,
                        change_entropy=-9.082223090817187e-7,
                        entropy_timederivative=-9.210638493017953e-7)
end

@testitem "linear_advection_FDSBP.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_strong_form.jl"),
                        l2=[0.00015746414238620223], linf=[0.0004289501170468002],
                        change_mass=2.7755575615628914e-17,
                        change_entropy=-9.082223090817187e-7,
                        entropy_timederivative=-9.210638493017953e-7)
end

@testitem "linear_advection_FDSBP_SAT.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl" with `surface_flux = (flux_central, flux_godunov)`
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_FDSBP_SAT.jl"),
                        boundary_conditions=boundary_condition_periodic,
                        l2=[0.00042725374698477996], linf=[0.0008162855376108169],
                        change_mass=-3.469446951953615e-16,
                        change_entropy=-4.298594954432744e-7,
                        entropy_timederivative=-6.420623096525803e-7)

    # Same errors as in "linear_advection.jl"
    D = couple_discontinuously(D_leg, uniform_mesh, Val(:minus))
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_FDSBP_SAT.jl"),
                        boundary_conditions=boundary_condition_periodic,
                        D=D,
                        l2=[0.000157464142385992], linf=[0.0004289501170479104],
                        change_mass=-2.7755575615628914e-17,
                        change_entropy=-9.082223083600738e-7,
                        entropy_timederivative=-9.21063849912418e-7)
end

@testitem "linear_advection_FD.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_FD.jl"),
                        l2=[0.03199259748601379], linf=[0.031230470786942688],
                        change_mass=0.00022039403436682758,
                        change_entropy=0.0038255363796113606,
                        entropy_timederivative=-0.0011976377220896395)

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
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_FD.jl"),
                            semi=semi,
                            l2=[0.00015746414238620223], linf=[0.0004289501170468002],
                            change_mass=2.7755575615628914e-17,
                            change_entropy=-9.082223090817187e-7,
                            entropy_timederivative=-9.210638493017953e-7)
    end
end

@testitem "linear_advection_flux_differencing.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                 "linear_advection_flux_differencing.jl"),
                        l2=[0.00015746414238611707], linf=[0.0004289501170472443],
                        change_mass=4.2327252813834093e-16,
                        change_entropy=-9.08222309914386e-7,
                        entropy_timederivative=-9.21063850300996e-7)
end

@testitem "linear_advection_flux_differencing_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection.jl"
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                 "linear_advection_flux_differencing_strong_form.jl"),
                        l2=[0.00015746414238603103], linf=[0.00042895011704957575],
                        change_mass=5.412337245047638e-16,
                        change_entropy=-9.082223089151853e-7,
                        entropy_timederivative=-9.210638502454849e-7)
end

@testitem "linear_advection_cfl.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                 "linear_advection_cfl.jl"),
                        l2=[0.00029528077087911356], linf=[0.0006703908593475028],
                        change_mass=-4.163336342344337e-17,
                        change_entropy=-0.0002495788935772403,
                        entropy_timederivative=-9.203964659865171e-7)

    @testset "linear_advection_cfl with time-varying CFL" begin
        # Same errors as above
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_cfl.jl"),
                            cfl=t -> 1.0,
                            l2=[0.00029528077087911356], linf=[0.0006703908593475028],
                            change_mass=-4.163336342344337e-17,
                            change_entropy=-0.0002495788935772403,
                            entropy_timederivative=-9.203964659865171e-7)
    end

    @testset "throw error if adaptive = true" begin
        @test_throws ArgumentError solve(ode, SSPRK53(), adaptive = true,
                                         callback = callbacks)
    end
end

@testitem "linear_advection_inhomogeneous_mesh.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                 "linear_advection_inhomogeneous_mesh.jl"),
                        l2=[0.00015505641092054123], linf=[0.0004295567744543316],
                        cons_error=[7.632783294297951e-17],
                        change_mass=7.632783294297951e-17,
                        change_entropy=-8.259112555530912e-7,
                        entropy_timederivative=-8.924542121424572e-7)

    coordinates = collect(coordinates_min:dx:coordinates_max)
    mesh = InhomogeneousMesh(coordinates)
    # Same errors as in "linear_advection.jl"
    @testset "Inhomogeneous mesh with homogeneous nodes" begin
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_inhomogeneous_mesh.jl"),
                            mesh=mesh,
                            l2=[0.0001574641423857773], linf=[0.00042895011704657815],
                            cons_error=[3.95516952522712e-16],
                            change_mass=3.95516952522712e-16,
                            change_entropy=-9.082223094147857e-7,
                            entropy_timederivative=-9.210638503565072e-7)
    end
end

@testitem "linear_advection_per_element.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION, "linear_advection_per_element.jl"),
                        l2=[0.0033676357851732484], linf=[0.009084487338805292],
                        cons_error=[7.494005416219807e-16],
                        change_mass=-7.494005416219807e-16,
                        change_entropy=-0.00025256281830143834,
                        entropy_timederivative=-0.0002625762280988875)

    Ds = [isodd(element) ? D_polydeg_2 : D_polydeg_4 for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                             volume_integral = VolumeIntegralStrongForm())
    @testset "Strong form" begin
        # Same errors as above
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_per_element.jl"),
                            solver=solver,
                            l2=[0.0033676357851732484], linf=[0.009084487338805292],
                            cons_error=[5.204170427930421e-16],
                            change_mass=5.204170427930421e-16,
                            change_entropy=-0.00025256281830349225,
                            entropy_timederivative=-0.0002625762280998867)
    end

    Ds = [isodd(element) ? D_polydeg_2 : D_polydeg_4 for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralWeakForm(flux_godunov),
                             volume_integral = VolumeIntegralFluxDifferencing())
    @testset "Flux differencing" begin
        # Same errors as above
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_per_element.jl"),
                            solver=solver,
                            l2=[0.0033676357851732484], linf=[0.009084487338805292],
                            cons_error=[3.95516952522712e-16],
                            change_mass=-3.95516952522712e-16,
                            change_entropy=-0.00025256281830293714,
                            entropy_timederivative=-0.0002625762281018851)
    end

    Ds = [isodd(element) ? D_polydeg_2 : D_polydeg_4 for element in eachelement(mesh)]
    solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                             volume_integral = VolumeIntegralFluxDifferencingStrongForm())
    @testset "Flux differencing, strong form" begin
        # Same errors as above
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_per_element.jl"),
                            solver=solver,
                            l2=[0.0033676357851732484], linf=[0.009084487338805292],
                            cons_error=[8.604228440844963e-16],
                            change_mass=8.604228440844963e-16,
                            change_entropy=-0.0002525628183037143,
                            entropy_timederivative=-0.0002625762281003863)
    end

    Ds = [legendre_derivative_operator(-1.0, 1.0, 4) for element in eachelement(mesh)]
    @testset "Same operators on each element" begin
        # Same errors as in "linear_advection.jl"
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_per_element.jl"),
                            Ds=Ds,
                            l2=[0.0001574641423857773], linf=[0.00042895011704657815],
                            cons_error=[6.453171330633722e-16],
                            change_mass=-6.453171330633722e-16,
                            change_entropy=-9.082223085821184e-7,
                            entropy_timederivative=-9.210638496903734e-7)
    end

    Ds = [legendre_derivative_operator(-1.0, 1.0, 4) for element in eachelement(mesh)]
    dx = element_spacing(mesh, 1)
    mesh = InhomogeneousMesh(collect(coordinates_min:dx:coordinates_max))
    @testset "Same operators on each element and InhomogeneousMesh" begin
        # Same errors as in "linear_advection.jl"
        @test_trixi_include(joinpath(EXAMPLES_DIR_ADVECTION,
                                     "linear_advection_per_element.jl"),
                            Ds=Ds, mesh=mesh,
                            l2=[0.0001574641423857773], linf=[0.00042895011704657815],
                            cons_error=[2.706168622523819e-16],
                            change_mass=-2.706168622523819e-16,
                            change_entropy=-9.082223089706964e-7,
                            entropy_timederivative=-9.210638503565072e-7)
    end
end
