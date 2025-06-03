using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 2.0
equations = LinearAdvectionEquation1D(advection_velocity)

function initial_condition_Gaussian(x, t, equations::LinearAdvectionEquation1D)
    x_trans = x - equations.advection_velocity * t
    return SVector(exp(-(x_trans + 0.5)^2 / 0.1))
end

# Create DG solver with polynomial degree = 3 and Godunov flux as surface flux
solver = DGSEM(polydeg = 3, surface_flux = flux_godunov)

coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0 # maximum coordinate

N_elements = 10 # number of elements
mesh = Mesh(coordinates_min, coordinates_max, N_elements)

boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition_Gaussian),
                       x_pos = boundary_condition_do_nothing)
# A semidiscretization collects data structures and functions for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition_Gaussian, solver;
                          boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_integrals = (mass, entropy))
callbacks = CallbackSet(analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            save_everystep = false, callback = callbacks, saveat = saveat)
