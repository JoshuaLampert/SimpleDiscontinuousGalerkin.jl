using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqSSPRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

# Create DG solver with polynomial degree = 3 and Godunov flux as surface flux
D = legendre_derivative_operator(-1.0, 1.0, 4)
solver = FDSBP(D, surface_flux = flux_godunov)

coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0 # maximum coordinate

N_elements = 10 # number of elements
mesh = Mesh(coordinates_min, coordinates_max, N_elements)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,))
stepsize_callback = StepsizeCallback(cfl = 1.0)
callbacks = CallbackSet(analysis_callback, summary_callback, stepsize_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, SSPRK53(), adaptive = false, dt = 1.0, # Will be overwritten by the `stepsize_callback`
            save_everystep = false, callback = callbacks, saveat = saveat)
