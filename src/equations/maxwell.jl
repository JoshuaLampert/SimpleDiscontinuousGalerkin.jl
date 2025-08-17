@doc raw"""
    MaxwellEquations1D(c = 299_792_458.0)

The Maxwell equations of electro dynamics
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
E \\ B
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
c^2 B \\ E
\end{pmatrix}
=
\begin{pmatrix}
0 \\ 0
\end{pmatrix}
```
in one dimension with speed of light `c = 299792458 m/s` (in vacuum).
In one dimension the Maxwell equations reduce to a wave equation.
The orthogonal magnetic (e.g.`B_y`) and electric field (`E_z`) propagate as waves
through the domain in `x`-direction.
For reference, see
- e.g. p.15 of Numerical Methods for Conservation Laws: From Analysis to Algorithms
  https://doi.org/10.1137/1.9781611975109

- or equation (1) in https://inria.hal.science/hal-01720293/document
"""
struct MaxwellEquations1D{RealT <: Real} <: AbstractEquations{1, 2}
    speed_of_light::RealT # c

    function MaxwellEquations1D(c::Real = 299_792_458.0)
        new{typeof(c)}(c)
    end
end

varnames(::typeof(cons2cons), ::MaxwellEquations1D) = ("E", "B")
varnames(::typeof(cons2prim), ::MaxwellEquations1D) = ("E", "B")

"""
    initial_condition_convergence_test(x, t, equations::MaxwellEquations1D)

A smooth initial condition used for convergence tests.
"""
function initial_condition_convergence_test(x, t, equations::MaxwellEquations1D)
    c = equations.speed_of_light
    char_pos = x - c * t
    char_min = x + c * t

    w1 = sinpi(2 * char_pos)
    w2 = sinpi(2 * char_min)

    E = (w1 + w2) / 2
    B = (w1 - w2) / (2 * c)
    return SVector(E, B)
end

# Calculate 1D flux for a single point
@inline function flux(u, equations::MaxwellEquations1D)
    E, B = u
    return SVector(equations.speed_of_light^2 * B, E)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed(u_ll, u_rr, equations::MaxwellEquations1D)
    return equations.speed_of_light
end

@inline cons2prim(u, ::MaxwellEquations1D) = u
@inline prim2cons(q, ::MaxwellEquations1D) = q

@inline function entropy(u, equations::MaxwellEquations1D)
    E, B = u
    return 0.5f0 * (E^2 + equations.speed_of_light^2 * B^2)
end
# Convert conservative variables to entropy variables
@inline function cons2entropy(u, equations::MaxwellEquations1D)
    E, B = u
    return SVector(E, equations.speed_of_light^2 * B)
end

function electric_field(u, equations::MaxwellEquations1D)
    return first(u)
end
function magnetic_field(u, equations::MaxwellEquations1D)
    return last(u)
end

pretty_form_utf(::typeof(electric_field)) = "∫E"
pretty_form_utf(::typeof(magnetic_field)) = "∫B"

function default_analysis_integrals(::MaxwellEquations1D)
    return (electric_field, magnetic_field, entropy, entropy_timederivative)
end

function (riemann_solver::RiemannSolver{MaxwellEquations1D{RealT}})(prob::RiemannProblem,
                                                                    xi) where {RealT}
    c = riemann_solver.equations.speed_of_light

    if xi < -c
        return prob.u_ll
    elseif xi > c
        return prob.u_rr
    else
        E_L, B_L = prob.u_ll
        E_R, B_R = prob.u_rr
        E_L_star = 0.5f0 * ((E_L + E_R) - c * (B_R - B_L))
        B_L_star = 0.5f0 * ((B_L + B_R) - (E_R - E_L) / c)
        return SVector(E_L_star, B_L_star)
    end
end
