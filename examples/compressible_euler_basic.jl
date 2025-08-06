using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations1D(1.4)

initial_condition = initial_condition_convergence_test
source_terms = source_terms_convergence_test

# Create DG solver with polynomial degree = 3 and Lax-Friedrichs flux as surface flux
surface_flux = flux_lax_friedrichs
volume_flux = flux_ranocha
solver = DGSEM(polydeg = 3, surface_integral = SurfaceIntegralWeakForm(surface_flux),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0 # maximum coordinate

N_elements = 10 # number of elements
mesh = Mesh(coordinates_min, coordinates_max, N_elements)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          source_terms = source_terms)

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
