using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK
using SummationByPartsOperatorsExtra: polynomialbases_derivative_operator

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0 # maximum coordinate

N_elements = 10 # number of elements
mesh = Mesh(coordinates_min, coordinates_max, N_elements)

# Create DG solver with alternating Gauss-Legendre and Gauss-Lobatto-Legendre operator
# with polynomial degree = 3 and Godunov flux as surface flux
D_GL = polynomialbases_derivative_operator(GaussLegendre, -1.0, 1.0, 4)
D_GLL = polynomialbases_derivative_operator(LobattoLegendre, -1.0, 1.0, 4)
Ds = [isodd(element) ? D_GL : D_GLL for element in eachelement(mesh)]
solver = PerElementFDSBP(Ds, surface_integral = SurfaceIntegralStrongForm(flux_godunov),
                         volume_integral = VolumeIntegralStrongForm())

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
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            save_everystep = false, callback = callbacks, saveat = saveat)
