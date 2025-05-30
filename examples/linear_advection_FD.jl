using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

coordinates_min = -2.0 # minimum coordinate
coordinates_max = 2.0 # maximum coordinate

D = derivative_operator(MattssonNordström2004(), 1, 2, -0.2, 0.2, 5)
solver = FDSBP(D, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
               volume_integral = VolumeIntegralStrongForm())

mesh = Mesh(coordinates_min, coordinates_max, 20)

# A semidiscretization collects data structures and functions for the spatial discretization
boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition),
                       x_pos = boundary_condition_do_nothing)
semi = Semidiscretization(mesh, equations, initial_condition, solver; boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
callbacks = CallbackSet(summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            save_everystep = false, callback = callbacks, saveat = saveat)
