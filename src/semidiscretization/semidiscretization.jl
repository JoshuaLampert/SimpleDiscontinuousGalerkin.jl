
"""
    Semidiscretization

A `struct` containing everything needed to describe a spatial semidiscretization
of an equation.
"""
struct Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                          SourceTerms, Solver, Cache}
    mesh::Mesh
    equations::Equations

    # This guy is a bit messy since we abuse it as some kind of "exact solution"
    # although this doesn't really exist...
    initial_condition::InitialCondition
    boundary_conditions::BoundaryConditions
    source_terms::SourceTerms
    solver::Solver
    cache::Cache
    function Semidiscretization{Mesh, Equations, InitialCondition, BoundaryConditions,
                                SourceTerms, Solver, Cache}(mesh::Mesh,
                                                            equations::Equations,
                                                            initial_condition::InitialCondition,
                                                            boundary_conditions::BoundaryConditions,
                                                            source_terms::SourceTerms,
                                                            solver::Solver,
                                                            cache::Cache) where {Mesh,
                                                                                 Equations,
                                                                                 InitialCondition,
                                                                                 BoundaryConditions,
                                                                                 SourceTerms,
                                                                                 Solver,
                                                                                 Cache}
        @assert ndims(mesh) == ndims(equations)

        new(mesh, equations, initial_condition, boundary_conditions, source_terms, solver,
            cache)
    end
end

"""
    Semidiscretization(mesh, equations, initial_condition, solver;
                       source_terms = nothing,
                       boundary_conditions = boundary_condition_periodic)

Construct a semidiscretization of a PDE.
"""
function Semidiscretization(mesh, equations, initial_condition, solver;
                            source_terms = nothing,
                            boundary_conditions = boundary_condition_periodic)
    _boundary_conditions = digest_boundary_conditions(boundary_conditions)
    _solver = digest_solver(mesh, solver)
    cache = (; create_cache(mesh, equations, _solver)...)
    Semidiscretization{typeof(mesh), typeof(equations), typeof(initial_condition),
                       typeof(_boundary_conditions), typeof(source_terms), typeof(_solver),
                       typeof(cache)}(mesh, equations, initial_condition,
                                      _boundary_conditions, source_terms,
                                      _solver, cache)
end

function Base.show(io::IO, semi::Semidiscretization)
    @nospecialize semi # reduce precompilation time

    print(io, "Semidiscretization(")
    print(io, semi.mesh)
    print(io, ", ", semi.equations)
    print(io, ", ", semi.initial_condition)
    print(io, ", ", semi.boundary_conditions)
    print(io, ", ", semi.source_terms)
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
        println(io, "    boundary condition: ", semi.boundary_conditions)
        print(io, "    source terms: ", semi.source_terms)
    end
end

@inline Base.ndims(semi::Semidiscretization) = ndims(semi.mesh)
@inline nvariables(semi::Semidiscretization) = nvariables(semi.equations)
@inline eachvariable(semi::Semidiscretization) = eachvariable(semi.equations)
@inline nelements(semi::Semidiscretization) = nelements(semi.mesh)
@inline eachelement(semi::Semidiscretization) = eachelement(semi.mesh)
@inline nnodes(semi::Semidiscretization, element) = nnodes(semi.solver, element)
@inline eachnode(semi::Semidiscretization, element) = eachnode(semi.solver, element)
@inline ndofs(semi::Semidiscretization) = ndofs(semi.mesh, semi.solver)
@inline Base.real(semi::Semidiscretization) = real(semi.mesh)

get_tmp_cache_scalar(cache) = cache.tmp_scalar

"""
    grid(semi)

Get the grid of a semidiscretization.
"""
grid(semi::Semidiscretization) = semi.cache.node_coordinates

"""
    flat_grid(semi)

Return a vector of the coordinates of all nodes in `semi`, flattened across all elements.
This is useful for plotting or other operations that require a single vector of coordinates.
"""
function flat_grid(semi::Semidiscretization)
    return collect(Iterators.flatten(parent(grid(semi))))
end

"""
    get_variable(u, v, semi)

Return the solution belonging to the variable `v` of the solution `u`
at one time step as a vector at every node across all elements.
"""
function get_variable(u, v, semi::Semidiscretization)
    get_variable(u, v, semi.solver)
end

# Here, `func` is a function that takes a vector at one element
# `u` is a vector of coefficients at all nodes of the element.
function integrate_on_element(func, u, D, element, jacobian)
    return jacobian[element] * integrate(func, u, D)
end
# This method is for integrating a vector quantity for all variables over the entire domain,
# such as the whole solution vector `u` (`Array{T, 3}` for DG methods with same basis across elements
# and `VectorOfArray{T, 3, Vector{Matrix{T}}}` for `PerElementFDSBP`).
function PolynomialBases.integrate(func,
                                   u::Union{Array{T, 3},
                                            VectorOfArray{T, 3, Vector{Matrix{T}}}},
                                   semi::Semidiscretization) where {T}
    integrals = zeros(real(semi), nvariables(semi))
    for v in eachvariable(semi)
        u_v = VectorOfArray([u[v, :, element] for element in eachelement(semi)])
        integrals[v] = integrate(func, u_v, semi)
    end
    return integrals
end

# This method is for integrating a scalar quantity over the entire domain.
function PolynomialBases.integrate(func, u, semi::Semidiscretization)
    return sum(integrate_on_element(func, u.u[element], get_basis(semi.solver, element),
                                    element, semi.cache.jacobian)
               for element in eachelement(semi))
end

function PolynomialBases.integrate(u, semi::Semidiscretization)
    integrate(identity, u, semi)
end

# Here, `func` is a function that takes a vector at one node and the equations, e.g., `mass` or `entropy`.
function integrate_quantity(func, u, semi::Semidiscretization)
    quantity = get_tmp_cache_scalar(semi.cache)
    integrate_quantity!(quantity, func, u, semi)
end

function compute_quantity!(quantity, func, u, mesh, equations, solver)
    for element in eachelement(mesh)
        for i in eachnode(solver, element)
            u_node = get_node_vars(u, equations, i, element)
            quantity[i, element] = func(u_node, equations)
        end
    end
end

function integrate_quantity!(quantity, func, u, semi::Semidiscretization)
    mesh, equations, solver, _ = mesh_equations_solver_cache(semi)
    compute_quantity!(quantity, func, u, mesh, equations, solver)
    integrate(quantity, semi)
end

function compute_quantity_timederivative!(quantity, func, du, u, mesh, equations, solver)
    for element in eachelement(mesh)
        for i in eachnode(solver, element)
            u_node = get_node_vars(u, equations, i, element)
            du_node = get_node_vars(du, equations, i, element)
            quantity[i, element] = dot(func(u_node, equations), du_node)
        end
    end
end

@inline function mesh_equations_solver_cache(semi::Semidiscretization)
    @unpack mesh, equations, solver, cache = semi
    return mesh, equations, solver, cache
end

function calc_error_norms(u, t, semi::Semidiscretization)
    calc_error_norms(u, t, semi.initial_condition, mesh_equations_solver_cache(semi)...)
end

function rhs!(du, u, semi::Semidiscretization, t)
    @unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

    @trixi_timeit timer() "rhs!" rhs!(du, u, t, mesh, equations, initial_condition,
                                      boundary_conditions, source_terms, solver, cache)

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

include("semidiscretization_overset_grid.jl")
