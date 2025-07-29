using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

initial_condition = initial_condition_convergence_test

# Create DG solver with polynomial degree = 3 and Godunov flux as surface flux
solver = DGSEM(polydeg = 3, surface_flux = flux_godunov)

a = -1.0
b = -0.1
c = 0.1
d = 1.0

N_elements = 10 # number of elements
mesh_left = Mesh(a, c, N_elements)
mesh_right = Mesh(b, d, N_elements)
mesh = OversetGridMesh(mesh_left, mesh_right)

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
stepsize_callback = StepsizeCallback(cfl = 0.5)
callbacks = CallbackSet(analysis_callback, summary_callback, stepsize_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, RDPK3SpFSAL49(), adaptive = false, dt = 1.0,
            save_everystep = false, callback = callbacks, saveat = saveat)

errs = errors(analysis_callback)
ints = integrals(analysis_callback)
using InteractiveUtils
clipboard("l2=$(errs.l2_error[:, end])\n,linf=$(errs.linf_error[:, end])\n,cons_error=$(errs.conservation_error[:, end])\n,change_mass=$(ints.mass[end] - ints.mass[1]),\nchange_entropy=$(ints.entropy[end] - ints.entropy[1]),\nentropy_timederivative=$(ints.entropy_timederivative[end])")