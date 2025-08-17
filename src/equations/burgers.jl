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
varnames(::typeof(cons2prim), ::BurgersEquation1D) = ("scalar",)

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
    omega = 2 * convert(RealT, pi) * f
    scalar = c + A * sin(omega * (x - t))

    return SVector(scalar)
end

"""
    source_terms_convergence_test(u, x, t, equations::BurgersEquation1D)

Source terms used for convergence tests in combination with
[`initial_condition_convergence_test`](@ref).
"""
@inline function source_terms_convergence_test(u, x, t, equations::BurgersEquation1D)
    # Same settings as in `initial_condition`
    RealT = eltype(x)
    c = 2
    A = 1
    L = 1
    f = 1.0f0 / L
    omega = 2 * convert(RealT, pi) * f
    du = omega * A * cos(omega * (x - t)) * (c - 1 + A * sin(omega * (x - t)))

    return SVector(du)
end

# Calculate 1D flux for a single point
@inline function flux(u, equation::BurgersEquation1D)
    return SVector(0.5f0 * u[1]^2)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed(u_ll, u_rr, equations::BurgersEquation1D)
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

function (riemann_solver::RiemannSolver{BurgersEquation1D})(prob::RiemannProblem, xi)
    u_ll, u_rr = prob.u_ll, prob.u_rr

    if u_ll[1] >= u_rr[1] # shock
        s = 0.5f0 * (u_ll[1] + u_rr[1])
        return xi < s ? u_ll : u_rr
    else # rarefaction wave
        if xi <= u_ll[1]
            return u_ll
        elseif xi >= u_rr[1]
            return u_rr
        else
            return xi
        end
    end
end
