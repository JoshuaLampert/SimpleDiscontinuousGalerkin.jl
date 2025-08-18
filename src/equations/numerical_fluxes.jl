"""
    flux_central(u_ll, u_rr, equations::AbstractEquations)

The classical central numerical flux `f((u_ll) + f(u_rr)) / 2`. When this flux is
used as volume flux, the discretization is equivalent to the classical weak form
DG method (except floating point errors).
"""
@inline function flux_central(u_ll, u_rr, equations::AbstractEquations)
    # Calculate regular 1D fluxes
    f_ll = flux(u_ll, equations)
    f_rr = flux(u_rr, equations)

    # Average regular fluxes
    return 0.5f0 * (f_ll + f_rr)
end

"""
    FluxPlusDissipation(numerical_flux, dissipation)

Combine a `numerical_flux` with a `dissipation` operator to create a new numerical flux.
"""
struct FluxPlusDissipation{NumericalFlux, Dissipation}
    numerical_flux::NumericalFlux
    dissipation::Dissipation
end

@inline function (numflux::FluxPlusDissipation)(u_ll, u_rr, equations)
    @unpack numerical_flux, dissipation = numflux

    return (numerical_flux(u_ll, u_rr, equations)
            +
            dissipation(u_ll, u_rr, equations))
end

function Base.show(io::IO, f::FluxPlusDissipation)
    print(io, "FluxPlusDissipation(", f.numerical_flux, ", ", f.dissipation, ")")
end

"""
    DissipationLocalLaxFriedrichs(max_abs_speed=max_abs_speed)

Create a local Lax-Friedrichs dissipation operator where the maximum absolute wave speed
is estimated as `max_abs_speed(u_ll, u_rr, equations)`, defaulting to [`max_abs_speed`](@ref).
"""
struct DissipationLocalLaxFriedrichs{MaxAbsSpeed}
    max_abs_speed::MaxAbsSpeed
end

DissipationLocalLaxFriedrichs() = DissipationLocalLaxFriedrichs(max_abs_speed)

@inline function (dissipation::DissipationLocalLaxFriedrichs)(u_ll, u_rr, equations)
    λ = dissipation.max_abs_speed(u_ll, u_rr, equations)
    return -0.5f0 * λ * (u_rr - u_ll)
end

function Base.show(io::IO, d::DissipationLocalLaxFriedrichs)
    print(io, "DissipationLocalLaxFriedrichs(", d.max_abs_speed, ")")
end

"""
    max_abs_speed(u_ll, u_rr, equations)

Simple and fast estimate of the maximal wave speed of the Riemann problem with left and right states
`u_ll, u_rr`, based only on the local wave speeds associated to `u_ll` and `u_rr`.
"""
function max_abs_speed end

max_abs_speeds(u_node, equations) = max_abs_speed(u_node, u_node, equations)

const FluxLaxFriedrichs{MaxAbsSpeed} = FluxPlusDissipation{typeof(flux_central),
                                                           DissipationLocalLaxFriedrichs{MaxAbsSpeed}}
"""
    FluxLaxFriedrichs(max_abs_speed=max_abs_speed)

Local Lax-Friedrichs (Rusanov) flux with maximum wave speed estimate provided by
`max_abs_speed`, cf. [`DissipationLocalLaxFriedrichs`](@ref) and
[`max_abs_speed`](@ref).
"""
function FluxLaxFriedrichs(max_abs_speed = max_abs_speed)
    FluxPlusDissipation(flux_central, DissipationLocalLaxFriedrichs(max_abs_speed))
end

function Base.show(io::IO, f::FluxLaxFriedrichs)
    print(io, "FluxLaxFriedrichs(", f.dissipation.max_abs_speed, ")")
end

"""
    flux_lax_friedrichs

See [`FluxLaxFriedrichs`](@ref).
"""
const flux_lax_friedrichs = FluxLaxFriedrichs()

"""
    FluxHLL(min_max_speed=min_max_speed_davis)

Create an HLL (Harten, Lax, van Leer) numerical flux where the minimum and maximum
wave speeds are estimated as
`λ_min, λ_max = min_max_speed(u_ll, u_rr, orientation_or_normal_direction, equations)`,
defaulting to [`min_max_speed_davis`](@ref).
Original paper:
- Amiram Harten, Peter D. Lax, Bram van Leer (1983)
  On Upstream Differencing and Godunov-Type Schemes for Hyperbolic Conservation Laws
  [DOI: 10.1137/1025002](https://doi.org/10.1137/1025002)
"""
struct FluxHLL{MinMaxSpeed}
    min_max_speed::MinMaxSpeed
end

FluxHLL() = FluxHLL(min_max_speed_davis)

"""
    min_max_speed_davis(u_ll, u_rr, equations)

Simple and fast estimates of the minimal and maximal wave speed of the Riemann problem with
left and right states `u_ll, u_rr`, usually based only on the local wave speeds associated to
`u_ll` and `u_rr`.

- S.F. Davis (1988)
  Simplified Second-Order Godunov-Type Methods
  [DOI: 10.1137/0909030](https://doi.org/10.1137/0909030)

See eq. (10.38) from
- Eleuterio F. Toro (2009)
  Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
  [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
See also [`FluxHLL`](@ref), [`min_max_speed_davis`](@ref).
"""
function min_max_speed_davis end

@inline function (numflux::FluxHLL)(u_ll, u_rr, equations)
    λ_min, λ_max = numflux.min_max_speed(u_ll, u_rr, equations)

    if λ_min >= 0 && λ_max >= 0
        return flux(u_ll, equations)
    elseif λ_max <= 0 && λ_min <= 0
        return flux(u_rr, equations)
    else
        f_ll = flux(u_ll, equations)
        f_rr = flux(u_rr, equations)
        inv_λ_max_minus_λ_min = inv(λ_max - λ_min)
        factor_ll = λ_max * inv_λ_max_minus_λ_min
        factor_rr = λ_min * inv_λ_max_minus_λ_min
        factor_diss = λ_min * λ_max * inv_λ_max_minus_λ_min
        return factor_ll * f_ll - factor_rr * f_rr + factor_diss * (u_rr - u_ll)
    end
end

Base.show(io::IO, numflux::FluxHLL) = print(io, "FluxHLL(", numflux.min_max_speed, ")")

"""
    flux_hll

See [`FluxHLL`](@ref).
"""
const flux_hll = FluxHLL()

"""
    RiemannProblem(u_ll, u_rr)

Create a Riemann problem with left and right states `u_ll` and `u_rr`,
which can be solved with a [`RiemannSolver`](@ref). This is used for the [`flux_godunov`](@ref)
numerical flux.
"""
struct RiemannProblem{ULType, URType}
    u_ll::ULType
    u_rr::URType
end

function Base.show(io::IO, prob::RiemannProblem)
    print(io, "RiemannProblem(", prob.u_ll, ", ", prob.u_rr, ")")
end

"""
    RiemannSolver(prob, equations)

An exact Riemann solver of the [`RiemannProblem`](@ref) prob for `equations`.
"""
struct RiemannSolver{Equations, ULType, URType, Cache}
    prob::RiemannProblem{ULType, URType}
    equations::Equations
    cache::Cache

    function RiemannSolver(prob::RiemannProblem{ULType, URType},
                           equations::AbstractEquations) where {ULType, URType}
        cache = nothing
        new{typeof(equations), ULType, URType, typeof(cache)}(prob, equations, cache)
    end
end

function Base.show(io::IO, riemann_solver::RiemannSolver)
    print(io, "RiemannSolver(", riemann_solver.prob, ", ", riemann_solver.equations, ")")
end

function (riemann_solver::RiemannSolver)(x, t)
    riemann_solver(x / t)
end

"""
    RiemannSolverSolution

A solution of a Riemann problem, which is a vector of states at different times.
Is returned when `solve`ing a [`RiemannProblem`](@ref) with a [`RiemannSolver`](@ref).
"""
struct RiemannSolverSolution{ULType, SolverType, XType, TType}
    solution::Vector{Vector{ULType}}
    solver::SolverType
    x::XType
    t::TType

    function RiemannSolverSolution(solution, riemann_solver::RiemannSolver, x, t)
        prob = riemann_solver.prob
        new{typeof(prob.u_ll), typeof(riemann_solver), typeof(x), typeof(t)}(solution,
                                                                             riemann_solver,
                                                                             x, t)
    end
end

Base.getindex(sol::RiemannSolverSolution, i::Int) = sol.solution[i]
function Base.setindex!(sol::RiemannSolverSolution, value, i::Int)
    setindex!(sol.solution, value, i)
end
Base.length(sol::RiemannSolverSolution) = length(sol.solution)

function Base.show(io::IO, sol::RiemannSolverSolution)
    print(io, "RiemannSolverSolution(", sol.solution, ", ", sol.solver, ", ", sol.x, ", ",
          sol.t, ")")
end

function Base.show(io::IO, ::MIME"text/plain", sol::RiemannSolverSolution)
    u_ll = sol.solver.prob.u_ll
    u_rr = sol.solver.prob.u_rr
    equations = sol.solver.equations
    print(io, "RiemannSolverSolution for ", equations, " with u_ll = ", u_ll, ", u_rr = ",
          u_rr, ", ", length(sol), " time steps, and ", length(sol.x), " spatial points")
end

# The additional kwarg `maxiters` is only used to make `@trixi_include` work, which tries to insert `maxiters` into `solve`
function CommonSolve.init(riemann_solver::RiemannSolver{Equations, ULType}, x, t;
                          maxiters = nothing) where {Equations, ULType}
    solution = Vector{Vector{ULType}}(undef, length(t))
    return RiemannSolverSolution(solution, riemann_solver, x, t)
end

function CommonSolve.solve!(sol::RiemannSolverSolution{ULType}) where {ULType}
    for (i, ti) in enumerate(sol.t)
        sol[i] = Vector{ULType}(undef, length(sol.x))
        for (j, xi) in enumerate(sol.x)
            sol[i][j] = sol.solver(xi, ti)
        end
    end
    return sol
end

"""
    flux_godunov(u_ll, u_rr, equations)

Numerical flux using the [`RiemannSolver`](@ref) to solve Riemann problems
with left and right states `u_ll` and `u_rr` for the given `equations` exactly.
"""
function flux_godunov(u_ll, u_rr, equations)
    prob = RiemannProblem(u_ll, u_rr)
    riemann_solver = RiemannSolver(prob, equations)
    godunov_state = riemann_solver(0)
    return flux(godunov_state, equations)
end
