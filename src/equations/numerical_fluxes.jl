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
    max_abs_speed_naive(u_ll, u_rr, equations)

Simple and fast estimate of the maximal wave speed of the Riemann problem with left and right states
`u_ll, u_rr`, based only on the local wave speeds associated to `u_ll` and `u_rr`.
"""
function max_abs_speed_naive end

"""
    max_abs_speed(u_ll, u_rr, equations)

Simple and fast estimate of the maximal wave speed of the Riemann problem with left and right states
`u_ll, u_rr`, based only on the local wave speeds associated to `u_ll` and `u_rr`.
Less diffusive, i.e., overestimating than [`max_abs_speed_naive`](@ref).
"""
@inline function max_abs_speed(u_ll, u_rr, equations)
    # Use naive version as "backup" if no specialized version for the equations at hand is available
    max_abs_speed_naive(u_ll, u_rr, equations)
end

max_abs_speeds(u_node, equations) = max_abs_speed(u_node, u_node, equations)

const FluxLaxFriedrichs{MaxAbsSpeed} = FluxPlusDissipation{typeof(flux_central),
                                                           DissipationLocalLaxFriedrichs{MaxAbsSpeed}}
"""
    FluxLaxFriedrichs(max_abs_speed=max_abs_speed)

Local Lax-Friedrichs (Rusanov) flux with maximum wave speed estimate provided by
`max_abs_speed`, cf. [`DissipationLocalLaxFriedrichs`](@ref) and
[`max_abs_speed_naive`](@ref).
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
    ExactRiemannSolver(equations)

An exact Riemann solver for the given `equations`.
"""
struct ExactRiemannSolver{Equations}
    u_ll::SVector
    u_rr::SVector
    equations::Equations

    function ExactRiemannSolver(u_ll::SVector{NVARS}, u_rr::SVector{NVARS},
                                equations::AbstractEquations{1, NVARS}) where {NVARS}
        new{typeof(equations)}(u_ll, u_rr, equations)
    end
end

function (riemann_solver::ExactRiemannSolver)(x, t)
    riemann_solver(x / t)
end

function flux_godunov(u_ll, u_rr, equations)
    riemann_solver = ExactRiemannSolver(u_ll, u_rr, equations)
    godunov_state = riemann_solver(zero(eltype(u_ll)))
    return flux(godunov_state, equations)
end
