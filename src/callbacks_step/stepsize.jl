# Mostly copied from Trixi.jl:
# https://github.com/trixi-framework/Trixi.jl/blob/e58392cae6683661579fc09a7398fd73a394381e/src/callbacks_step/stepsize.jl
"""
    StepsizeCallback(; cfl=1.0, interval = 1)

Set the time step size according to a CFL condition with CFL number `cfl`
if the time integration method isn't adaptive itself.

The supplied keyword argument `cfl` must be either a `Real` number or
a function of time `t` returning a `Real` number.
By default, the timestep will be adjusted at every step.
For different values of `interval`, the timestep will be adjusted every `interval` steps.
"""
mutable struct StepsizeCallback{CflType}
    cfl_number::CflType
    interval::Int
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:StepsizeCallback})
    @nospecialize cb # reduce precompilation time

    stepsize_callback = cb.affect!
    @unpack cfl_number, interval = stepsize_callback
    print(io, "StepsizeCallback(",
          "cfl_number=", cfl_number, ", ",
          "interval=", interval, ")")
end

function StepsizeCallback(; cfl = 1.0, interval = 1)
    stepsize_callback = StepsizeCallback{typeof(cfl)}(cfl, interval)

    DiscreteCallback(stepsize_callback, stepsize_callback, # the first one is the condition, the second the affect!
                     save_positions = (false, false),
                     initialize = initialize!)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u, t,
                     integrator) where {Condition, Affect! <: StepsizeCallback}
    cb.affect!(integrator)
end

# this method is called to determine whether the callback should be activated
function (stepsize_callback::StepsizeCallback)(u, t, integrator)
    @unpack interval = stepsize_callback

    # Although the CFL-based timestep is usually not used with
    # adaptive time integration methods, we still check the accepted steps `naccept` here.
    return interval > 0 && integrator.stats.naccept % interval == 0
end

# This method is called as callback during the time integration.
@inline function (stepsize_callback::StepsizeCallback)(integrator)
    if integrator.opts.adaptive
        throw(ArgumentError("The `StepsizeCallback` has no effect when using an adaptive time integration scheme. Please remove the `StepsizeCallback` or set `adaptive = false` in `solve`."))
    end

    t = integrator.t
    u = integrator.u
    semi = integrator.p
    @unpack cfl_number = stepsize_callback

    # Dispatch based on semidiscretization
    dt = calculate_dt(u, t, cfl_number, semi)

    set_proposed_dt!(integrator, dt)
    integrator.opts.dtmax = dt
    integrator.dtcache = dt

    # avoid re-evaluating possible FSAL stages
    u_modified!(integrator, false)
    return nothing
end

# General case for a single (i.e., non-coupled) semidiscretization
# Case for constant `cfl_number`.
function calculate_dt(u, t, cfl_number::Real, semi)
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)

    dt = cfl_number * max_dt(u, t, mesh, equations, solver, cache)
    return dt
end
# Case for `cfl_number` as a function of time `t`.
function calculate_dt(u, t, cfl_number, semi)
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)

    dt = cfl_number(t) * max_dt(u, t, mesh, equations, solver, cache)
    return dt
end

function max_dt(u, t, mesh, equations, dg, cache)
    # to avoid a division by zero if the speed vanishes everywhere,
    # e.g. for steady-state linear advection
    max_scaled_speed = nextfloat(zero(t))

    for element in eachelement(mesh)
        max_lambda1 = zero(max_scaled_speed)
        for i in eachnode(dg, element)
            u_node = get_node_vars(u, equations, i, element)
            lambda1, = max_abs_speeds(u_node, equations)
            max_lambda1 = max(max_lambda1, lambda1)
        end
        # We need to take the length of the basis into account because
        # in contrast to Trixi.jl, it is not always on the reference element [-1, 1].
        basis = get_basis(dg, element)
        dx_basis = last(grid(basis)) - first(grid(basis))
        inv_jacobian = 1 / (dx_basis * cache.jacobian[element])
        max_scaled_speed = max(max_scaled_speed, inv_jacobian * max_lambda1)
    end

    return 1 / (nnodes(dg) * max_scaled_speed)
end
