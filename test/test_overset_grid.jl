@testitem "linear_advection_overset_grid.jl" setup=[Setup] begin
    # Not mass conservative because we miss integrating the part from the left boundary of the left
    # overlap element to b.
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_overset_grid.jl"),
                        l2=[2.3825331813104963e-5], linf=[3.095999626845369e-5],
                        cons_error=[6.721897816142075e-9],
                        change_mass=-6.721897816142075e-9,
                        change_entropy=8.533662195886471e-6,
                        entropy_timederivative=0.09450660830795257)

    dx = (c - a) / N_elements
    coordinates = collect(a:dx:c)
    push!(coordinates, b)
    mesh_left = InhomogeneousMesh(coordinates)
    # Adding b to the left mesh yields mass conservation again
    @testset "InhomogeneousMesh" begin
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_overset_grid.jl"),
                            mesh_left=mesh_left,
                            l2=[2.0422480731399535e-5], linf=[3.951360818943428e-5],
                            cons_error=[2.220446049250313e-16],
                            change_mass=2.220446049250313e-16,
                            change_entropy=-1.9940662454587255e-8,
                            entropy_timederivative=-1.4171278622798766e-8)
    end

    function initial_condition_Gaussian(x, t, equations::LinearAdvectionEquation1D)
        x_trans = x - equations.advection_velocity * t
        return SVector(exp(-(x_trans + 0.5)^2 / 0.1))
    end
    initial_condition = initial_condition_Gaussian
    boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition),
                           x_pos = boundary_condition_do_nothing)
    semi = Semidiscretization(mesh, equations, initial_condition, solver;
                              boundary_conditions)
    @testset "Boundary conditions" begin
        @test_trixi_include(joinpath(examples_dir(), "linear_advection_overset_grid.jl"),
                            semi=semi,
                            l2=[9.287766266389129e-6], linf=[4.926556788237193e-5],
                            cons_error=[0.5462913980698885],
                            change_mass=-0.5462913980698885,
                            change_entropy=-0.19785618067708552,
                            entropy_timederivative=-0.006735287186087797)
    end
end

@testitem "linear_advection_overset_grid_strong_form.jl" setup=[Setup] begin
    # Same errors as in "linear_advection_overset_grid.jl"
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_overset_grid_strong_form.jl"),
                        l2=[2.3825331814026118e-5], linf=[3.095999626745449e-5],
                        cons_error=[6.721899148409705e-9],
                        change_mass=-6.721899148409705e-9,
                        change_entropy=8.53366219633056e-6,
                        entropy_timederivative=0.09450660830795288)

    @testset "VolumeIntegralFluxDifferencing" begin
        # Same errors as in "linear_advection_overset_grid.jl"
        @test_trixi_include(joinpath(examples_dir(),
                                     "linear_advection_overset_grid_strong_form.jl"),
                            surface_integral=SurfaceIntegralWeakForm(flux_godunov),
                            volume_integral=VolumeIntegralFluxDifferencing(flux_central),
                            l2=[2.382533181422071e-5], linf=[3.095999626634427e-5],
                            cons_error=[6.72189803818668e-9],
                            change_mass=-6.72189803818668e-9,
                            change_entropy=8.533662197107716e-6,
                            entropy_timederivative=0.0945066083079539)
    end

    @testset "VolumeIntegralFluxDifferencingStrongForm" begin
        # Same errors as in "linear_advection_overset_grid.jl"
        @test_trixi_include(joinpath(examples_dir(),
                                     "linear_advection_overset_grid_strong_form.jl"),
                            surface_integral=SurfaceIntegralStrongForm(flux_godunov),
                            volume_integral=VolumeIntegralFluxDifferencingStrongForm(flux_central),
                            l2=[2.3825331813872226e-5], linf=[3.095999626678836e-5],
                            cons_error=[6.72189803818668e-9],
                            change_mass=-6.72189803818668e-9,
                            change_entropy=8.533662195886471e-6,
                            entropy_timederivative=0.09450660830795123)
    end
end

@testitem "linear_advection_overset_grid_cfl.jl" setup=[Setup] begin
    # Not mass conservative because we miss integrating the part from the left boundary of the left
    # overlap element to b.
    @test_trixi_include(joinpath(examples_dir(), "linear_advection_overset_grid_cfl.jl"),
                        l2=[2.381442225978849e-5], linf=[3.098598080530923e-5],
                        cons_error=[2.1519649284762465e-8],
                        change_mass=-2.1519649284762465e-8,
                        change_entropy=8.54207263600859e-6,
                        entropy_timederivative=0.09450669932241434)
end

@testitem "linear_advection_overset_grid_per_element.jl" setup=[Setup] begin
    # Mass conservative because b is exactly in the left mesh (11 elements).
    @test_trixi_include(joinpath(examples_dir(),
                                 "linear_advection_overset_grid_per_element.jl"),
                        l2=[0.0007287829317534878], linf=[0.0016309648551183775],
                        cons_error=[4.440892098500626e-16],
                        change_mass=4.440892098500626e-16,
                        change_entropy=-1.049475386327714e-5,
                        entropy_timederivative=-1.0654816822358582e-5)
end

@testitem "Maxwell_overset_grid.jl" setup=[Setup] begin
    # Not mass conservative because we miss integrating the part from the left boundary of the left
    # overlap element to b.
    @test_trixi_include(joinpath(examples_dir(), "Maxwell_overset_grid.jl"),
                        l2=[0.00023060901614886724, 0.00026184365443377735],
                        linf=[0.0002498594703230772, 0.00036105557236677105],
                        cons_error=[0.001469440514095803, 0.024591060510140168],
                        change_entropy=0.0032048130928074015,
                        entropy_timederivative=0.12743862294720792)
end

@testitem "Burgers_overset_grid.jl" setup=[Setup] begin
    # Mass conservative because we choose 11 elements meaning b is exactly an interface.
    @test_trixi_include(joinpath(examples_dir(), "Burgers_overset_grid.jl"),
                        l2=[1.0898127073451138], linf=[0.7933574712125537],
                        cons_error=[6.159517340620368e-13],
                        change_mass=-6.159517340620368e-13,
                        change_entropy=-0.43769950084048936,
                        entropy_timederivative=-0.08354232871802136)
end
