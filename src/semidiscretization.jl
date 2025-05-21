
"""
    Semidiscretization

A `struct` containing everything needed to describe a spatial semidiscretization
of an equation.
"""
struct Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                          Solver, Cache}
    mesh::Mesh
    equations::Equations

    # This guy is a bit messy since we abuse it as some kind of "exact solution"
    # although this doesn't really exist...
    initial_condition::InitialCondition
    boundary_conditions::BoundaryConditions
    solver::Solver
    cache::Cache
    function Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                                Solver, Cache}(mesh::Mesh, equations::Equations,
                                               initial_condition::InitialCondition,
                                               boundary_conditions::BoundaryConditions,
                                               solver::Solver,
                                               cache::Cache) where {Mesh, Equations,
                                                                    InitialCondition,
                                                                    BoundaryConditions,
                                                                    Solver, Cache}
        @assert ndims(mesh) == ndims(equations)

        new(mesh, equations, initial_condition, boundary_conditions, solver, cache)
    end
end

"""
    Semidiscretization(mesh, equations, initial_condition, solver;
                       boundary_conditions = boundary_condition_periodic)

Construct a semidiscretization of a PDE.
"""
function Semidiscretization(mesh, equations, initial_condition, solver;
                            boundary_conditions = boundary_condition_periodic)
    cache = create_cache(mesh, equations, solver, initial_condition,
                         boundary_conditions)
    _boundary_conditions = digest_boundary_conditions(boundary_conditions, mesh, solver,
                                                      cache)
    Semidiscretization{typeof(mesh), typeof(equations), typeof(initial_condition),
                       typeof(_boundary_conditions), typeof(solver), typeof(cache)}(mesh,
                                                                                    equations,
                                                                                    initial_condition,
                                                                                    _boundary_conditions,
                                                                                    solver,
                                                                                    cache)
end

function Base.show(io::IO, semi::Semidiscretization)
    @nospecialize semi # reduce precompilation time

    print(io, "Semidiscretization(")
    print(io, semi.mesh)
    print(io, ", ", semi.equations)
    print(io, ", ", semi.initial_condition)
    print(io, ", ", semi.boundary_conditions)
    print(io, ", ", semi.solver)
    print(io, "))")
end

function Base.show(io::IO, ::MIME"text/plain", semi::Semidiscretization)
    @nospecialize semi # reduce precompilation time

    if get(io, :compact, false)
        show(io, semi)
    else
        println(io, "Semidiscretization")
        println(io, "    #spatial dimensions: ", ndims(semi))
        println(io, "    mesh: ", semi.mesh)
        println(io, "    equations: ", get_name(semi.equations))
        println(io, "    initial condition: ", semi.initial_condition)
        print(io, "    boundary condition: ", semi.boundary_conditions)
    end
end

@inline Base.ndims(semi::Semidiscretization) = ndims(semi.mesh)
@inline nvariables(semi::Semidiscretization) = nvariables(semi.equations)
@inline eachvariable(semi::Semidiscretization) = eachvariable(semi.equations)
@inline nelements(semi::Semidiscretization) = nelements(semi.mesh)
@inline eachelement(semi::Semidiscretization) = eachelement(semi.mesh)
@inline ndofs(semi::Semidiscretization) = ndofs(semi.mesh, semi.solver)
@inline Base.real(semi::Semidiscretization) = real(semi.solver)

"""
    grid(semi)

Get the grid of a semidiscretization.
"""
SummationByPartsOperators.grid(semi::Semidiscretization) = semi.cache.node_coordinates

@inline function mesh_equations_solver_cache(semi::Semidiscretization)
    @unpack mesh, equations, solver, cache = semi
    return mesh, equations, solver, cache
end

function rhs!(du, u, semi::Semidiscretization, t)
    @unpack mesh, equations, initial_condition, boundary_conditions, solver, cache = semi

    @trixi_timeit timer() "rhs!" rhs!(du, u, t, mesh, equations, initial_condition,
                                      boundary_conditions, solver, cache)

    return nothing
end

function compute_coefficients(func, t, semi::Semidiscretization)
    u = allocate_coefficients(mesh_equations_solver_cache(semi)...)
    compute_coefficients!(u, func, t, semi)
    return u
end

function compute_coefficients!(u, func, t, semi::Semidiscretization)
    compute_coefficients!(u, func, t, mesh_equations_solver_cache(semi)...)
end

"""
    semidiscretize(semi::Semidiscretization, tspan)

Wrap the semidiscretization `semi` as an ODE problem in the time interval `tspan`
that can be passed to `solve` from the [SciML ecosystem](https://diffeq.sciml.ai/latest/).
"""
function semidiscretize(semi::Semidiscretization, tspan)
    u0 = compute_coefficients(semi.initial_condition, first(tspan), semi)
    iip = true # is-inplace, i.e., we modify a vector when calling rhs!
    return ODEProblem{iip}(rhs!, u0, tspan, semi)
end
