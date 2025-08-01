@doc raw"""
    BurgersEquation1D

The inviscid Burgers' equation
```math
\partial_t u + \frac{1}{2} \partial_1 u^2 = 0
```
in one space dimension.
"""
struct BurgersEquation1D <: AbstractEquations{1, 1} end

varnames(::typeof(cons2cons), ::BurgersEquation1D) = ("scalar",)

"""
    initial_condition_convergence_test(x, t, equations::BurgersEquation1D)

A smooth initial condition used for convergence tests.
"""
function initial_condition_convergence_test(x, t, equation::BurgersEquation1D)
    RealT = eltype(x)
    c = 2
    A = 1
    L = 1
    f = 1.0f0 / L
    omega = 2 * f
    scalar = c + A * sinpi(omega * (x - t))

    return SVector(scalar)
end

# Calculate 1D flux for a single point
@inline function flux(u, equation::BurgersEquation1D)
    return SVector(0.5f0 * u[1]^2)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr,
                                     equations::BurgersEquation1D)
    u_L = u_ll[1]
    u_R = u_rr[1]

    return max(abs(u_L), abs(u_R))
end

@doc raw"""
    flux_ec(u_ll, u_rr, equations::BurgersEquation1D)

Entropy-conserving, symmetric flux for the inviscid Burgers' equation.
```math
F(u_L, u_R) = \frac{u_L^2 + u_L u_R + u_R^2}{6}
```
"""
function flux_ec(u_ll, u_rr, equation::BurgersEquation1D)
    u_L = u_ll[1]
    u_R = u_rr[1]

    return SVector((u_L^2 + u_L * u_R + u_R^2) / 6)
end

"""
    flux_godunov(u_ll, u_rr, equations::BurgersEquation1D)

Godunov (upwind) numerical flux for the inviscid Burgers' equation.
See https://metaphor.ethz.ch/x/2019/hs/401-4671-00L/literature/mishra_hyperbolic_pdes.pdf ,
section 4.1.5 and especially equation (4.16).
"""
function flux_godunov(u_ll, u_rr, equation::BurgersEquation1D)
    u_L = u_ll[1]
    u_R = u_rr[1]

    return SVector(0.5f0 * max(max(u_L, 0)^2, min(u_R, 0)^2))
end
