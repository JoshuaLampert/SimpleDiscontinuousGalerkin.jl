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

N_elements = 10 # number of elements
mesh_left = Mesh(a, c, N_elements)
mesh_right = Mesh(b, d, N_elements)
mesh = OversetGridMesh(mesh_left, mesh_right)

mesh_left = Mesh(a, c, 1)
mesh_right = Mesh(b, d, 40)
mesh = OversetGridMesh(mesh_left, mesh_right)

x_L_ref = -1.0
x_R_ref = 1.0
p = 4
# n_nodes = 48
N_elements = 51
n_nodes = N_elements
D_FD = derivative_operator(MattssonNordström2004(), 1, p, x_L_ref, x_R_ref, n_nodes)

surface_integral = SurfaceIntegralStrongForm(flux_godunov)
volume_integral = VolumeIntegralStrongForm()

Ds_left = [D_FD for element in eachelement(mesh_left)]
solver_left = PerElementFDSBP(Ds_left,
                              surface_integral = surface_integral,
                              volume_integral = volume_integral)

D_GLL = legendre_derivative_operator(x_L_ref, x_R_ref, 3 + 1)
Ds_right = [D_GLL for element in eachelement(mesh_right)]
solver_right = PerElementFDSBP(Ds_right,
                               surface_integral = surface_integral,
                               volume_integral = volume_integral)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, (solver_left, solver_right))

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; io = devnull, interval = 10,
                                     extra_analysis_errors = (:conservation_error,))
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            save_everystep = false, callback = callbacks, saveat = saveat)
