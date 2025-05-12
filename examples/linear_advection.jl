using SimpleDiscontinuousGalerkin
using OrdinaryDiffEqLowStorageRK

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = 1.0
equations = LinearAdvectionEquation1D(advection_velocity)

function initial_condition(x, t, equations::LinearAdvectionEquation1D)
    x_trans = x - equations.advection_velocity * t
    return SVector(sinpi(x_trans))
end

# Create DG solver with polynomial degree = 3 and (local) Godunov flux as surface flux
solver = DGSEM(polydeg = 3, surface_flux = flux_godunov)

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
saveat = range(tspan..., length = 100)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            saveat = saveat)
