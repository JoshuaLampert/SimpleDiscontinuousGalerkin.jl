using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0 # maximum coordinate

# We couple the discontinuous Legendre derivative operator on a uniform mesh
# directly with the tools from SummationByPartsOperators.jl to create a
# global DG operator and we only use one big element here. This is equivalent to
# using the DGSEM solver with polynomial degree 3 and a central flux for the interior
# and a Godunov flux for the boundary. This is equivalent to the global SBP-SAT method
# du = -a * (D * u) - aM^{-1} * e_L * (e_L^T * u - g(t)),
# where D is a discontinuously coupled Legendre derivative operator.
p = 3
D_leg = legendre_derivative_operator(-1.0, 1.0, p + 1)
N_elements = 10
uniform_mesh = UniformMesh1D(coordinates_min, coordinates_max, N_elements)
D = couple_discontinuously(D_leg, uniform_mesh)
# The interior flux doesn't matter here because we only have one element.
solver = FDSBP(D, surface_integral = SurfaceIntegralStrongForm(flux_central, flux_godunov),
               volume_integral = VolumeIntegralStrongForm())

mesh = Mesh(coordinates_min, coordinates_max, 1) # use only one element because we already have a global operator

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
