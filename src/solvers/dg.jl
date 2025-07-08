"""
    DG(; basis, surface_integral, volume_integral)

Create a discontinuous Galerkin method.
If `basis isa LegendreDerivativeOperator`, this creates a [`DGSEM`](@ref).
"""
struct DG{Basis, SurfaceIntegral, VolumeIntegral}
    basis::Basis
    surface_integral::SurfaceIntegral
    volume_integral::VolumeIntegral
end

function Base.show(io::IO, solver::DG)
    @nospecialize solver # reduce precompilation time

    print(io, "DG{", real(solver), "}(")
    print(io, solver.basis)
    print(io, ", ", solver.surface_integral)
    print(io, ", ", solver.volume_integral)
    print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", solver::DG)
    @nospecialize solver # reduce precompilation time

    if get(io, :compact, false)
        show(io, solver)
    else
        println(io, "DG{" * string(real(solver)) * "}")
        println(io, "    basis: ", solver.basis)
        println(io, "    surface integral: ", solver.surface_integral |> typeof |> nameof)
        print(io, "    volume integral: ", solver.volume_integral |> typeof |> nameof)
    end
end

Base.summary(io::IO, solver::DG) = print(io, "DG(" * summary(solver.basis) * ")")

@inline Base.real(solver::DG) = real(solver.basis)

# This method is only supported for solvers with a fixed number of nodes per element.
grid(solver::DG) = grid(solver.basis)
grid(solver::DG, element) = grid(solver.basis)

"""
    eachnode(solver::DG, element)

Return an iterator over the indices that specify the location in relevant data structures
for the nodes in a specific `element` in `solver`.
In particular, not the nodes themselves are returned.
"""
@inline eachnode(solver::DG, element) = Base.OneTo(nnodes(solver, element))
# This method is only supported for solvers with a fixed number of nodes per element.
@inline nnodes(solver::DG) = length(grid(solver))
@inline nnodes(solver::DG, element) = length(grid(solver, element))
@inline function ndofs(mesh::AbstractMesh, solver::DG)
    sum(nnodes(solver, element) for element in eachelement(mesh))
end

@inline function get_node_coords(x, equations, ::DG, indices...)
    return x[indices...]
end

# Adapted from Trixi.jl
# https://github.com/trixi-framework/Trixi.jl/blob/75d8c67629562efd24b2a04e46d22b0a1f4f572c/src/solvers/dg.jl#L539
@inline function get_node_vars(u, equations, indices...)
    # There is a cut-off at `n == 10` inside of the method
    # `ntuple(f::F, n::Integer) where F` in Base at ntuple.jl:17
    # in Julia `v1.5`, leading to type instabilities if
    # more than ten variables are used. That's why we use
    # `Val(...)` below.
    # We use `@inline` to make sure that the `getindex` calls are
    # really inlined, which might be the default choice of the Julia
    # compiler for standard `Array`s but not necessarily for more
    # advanced array types such as `PtrArray`s, cf.
    # https://github.com/JuliaSIMD/VectorizationBase.jl/issues/55
    SVector(ntuple(@inline(v->u[v, indices...]), Val(nvariables(equations))))
end

@inline function set_node_vars!(u, u_node, equations, indices...)
    for v in eachvariable(equations)
        u[v, indices...] = u_node[v]
    end
    return nothing
end

"""
    get_variable(u, v, ::DG)

Return the solution belonging to the variable `v` of the solution `u`
at one time step as a vector at every node across all elements.
"""
function get_variable(u, v, ::DG)
    return vec(u[v, :, :])
end

function allocate_coefficients(mesh::AbstractMesh, equations, solver::DG, cache)
    return allocate_coefficients(mesh, equations, solver)
end

function allocate_coefficients(mesh::AbstractMesh, equations, solver::DG)
    return zeros(real(solver), nvariables(equations), nnodes(solver), nelements(mesh))
end

function compute_coefficients!(u, func, t, mesh::AbstractMesh, equations, solver::DG, cache)
    compute_coefficients!(u, func, t, mesh, equations, solver, cache, cache.node_coordinates)
end

function compute_coefficients!(u, func, t, mesh::AbstractMesh, equations, solver::DG,
                               cache, node_coordinates)
    for element in eachelement(mesh)
        for i in eachnode(solver, element)
            x_node = get_node_coords(node_coordinates, equations, solver, i, element)
            u_node = func(x_node, t, equations)
            set_node_vars!(u, u_node, equations, i, element)
        end
    end
end

function reset_du!(du)
    du .= zero(du)
end

function rhs!(du, u, t, mesh::AbstractMesh, equations, initial_condition,
              boundary_conditions, solver::DG, cache)
    @trixi_timeit timer() "reset ∂u/∂t" reset_du!(du)

    @trixi_timeit timer() "volume integral" begin
        calc_volume_integral!(du, u, mesh, equations,
                              solver.volume_integral, solver, cache)
    end

    @trixi_timeit timer() "interface flux" begin
        calc_interface_flux!(cache.surface_flux_values, u, mesh,
                             equations, solver.surface_integral, solver, cache)
    end

    @trixi_timeit timer() "boundary flux" begin
        calc_boundary_flux!(cache.surface_flux_values, u, t, boundary_conditions, mesh,
                            equations, solver.surface_integral, solver)
    end

    @trixi_timeit timer() "surface integral" begin
        calc_surface_integral!(du, u, mesh, equations,
                               solver.surface_integral, solver, cache)
    end

    @trixi_timeit timer() "Jacobian" apply_jacobian!(du, mesh, equations, solver, cache)

    return nothing
end
