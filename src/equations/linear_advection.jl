@doc raw"""
    LinearAdvectionEquation1D

The linear scalar advection equation
```math
\partial_t u + a \partial_1 u  = 0
```
in one space dimension with constant velocity `a`.
"""
struct LinearAdvectionEquation1D{RealT <: Real} <: AbstractEquations{1, 1}
    advection_velocity::RealT
end

varnames(::typeof(cons2cons), ::LinearAdvectionEquation1D) = ("u",)
varnames(::typeof(cons2prim), ::LinearAdvectionEquation1D) = ("u",)

"""
    initial_condition_convergence_test(x, t, equations::LinearAdvectionEquation1D)

A smooth initial condition used for convergence tests.
"""
function initial_condition_convergence_test(x, t, equations::LinearAdvectionEquation1D)
    x_trans = x - equations.advection_velocity * t
    return SVector(sinpi(x_trans))
end

@inline function flux(u, equation::LinearAdvectionEquation1D)
    a = equation.advection_velocity
    return a * u
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed(u_ll, u_rr, equation::LinearAdvectionEquation1D)
    return abs(equation.advection_velocity)
end

function (riemann_solver::RiemannSolver{LinearAdvectionEquation1D{RealT}})(xi) where {RealT}
    u_ll, u_rr = riemann_solver.prob.u_ll, riemann_solver.prob.u_rr
    a = riemann_solver.equations.advection_velocity
    return a >= xi ? u_ll : u_rr
end
