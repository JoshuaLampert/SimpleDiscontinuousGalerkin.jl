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
@inline function max_abs_speed_naive(u_ll, u_rr,
                                     equation::LinearAdvectionEquation1D)
    return abs(equation.advection_velocity)
end

"""
    flux_godunov(u_ll, u_rr, equations::LinearAdvectionEquation1D)

Godunov (upwind) flux for the 1D linear scalar advection equation.
Essentially first order upwind, see e.g. https://math.stackexchange.com/a/4355076/805029 .
"""
function flux_godunov(u_ll, u_rr, equation::LinearAdvectionEquation1D)
    u_L = u_ll[1]
    u_R = u_rr[1]

    a = equation.advection_velocity
    if a >= 0
        return SVector(a * u_L)
    else
        return SVector(a * u_R)
    end
end
