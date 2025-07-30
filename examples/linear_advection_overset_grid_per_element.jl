using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

a = -1.0
b = -0.1
c = 0.1
d = 1.0

mesh_left = Mesh(a, c, 11)
mesh_right = Mesh(b, d, 10)
mesh = OversetGridMesh(mesh_left, mesh_right)

# Use alternating second-order and fourth-order Legendre derivative operators.
# Wrap them in `Ref` to allow sharing them between elements without copying.
D_polydeg_2 = legendre_derivative_operator(-1.0, 1.0, 3)
D_polydeg_4 = legendre_derivative_operator(-1.0, 1.0, 5)
Ds_left = [isodd(element) ? D_polydeg_4 : D_polydeg_2 for element in eachelement(mesh_left)]
solver_left = PerElementFDSBP(Ds_left,
                              surface_integral = SurfaceIntegralWeakForm(flux_godunov),
                              volume_integral = VolumeIntegralWeakForm())
Ds_right = [isodd(element) ? D_polydeg_4 : D_polydeg_2
            for element in eachelement(mesh_right)]
solver_right = PerElementFDSBP(Ds_right,
                               surface_integral = SurfaceIntegralWeakForm(flux_godunov),
                               volume_integral = VolumeIntegralWeakForm())

# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, (solver_left, solver_right))

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
