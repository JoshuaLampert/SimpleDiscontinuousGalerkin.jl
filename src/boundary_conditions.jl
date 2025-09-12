function digest_boundary_conditions(boundary_conditions)
    (; x_neg = boundary_conditions, x_pos = boundary_conditions)
end

function digest_boundary_conditions(boundary_conditions::NTuple{2, Any})
    (; x_neg = boundary_conditions[1], x_pos = boundary_conditions[2])
end
function digest_boundary_conditions(boundary_conditions::NamedTuple{Keys, ValueTypes}) where {
                                                                                              Keys,
                                                                                              ValueTypes <:
                                                                                              NTuple{2,
                                                                                                     Any}
                                                                                              }
    @unpack x_neg, x_pos = boundary_conditions
    (; x_neg, x_pos)
end

struct BoundaryConditionPeriodic end

"""
    boundary_condition_periodic = BoundaryConditionPeriodic()

A singleton struct indicating periodic boundary conditions.
"""
const boundary_condition_periodic = BoundaryConditionPeriodic()

function Base.show(io::IO, ::BoundaryConditionPeriodic)
    print(io, "boundary_condition_periodic")
end

@inline function (::BoundaryConditionPeriodic)(u, x, t, mesh, equations, solver, is_left, cache)
    N_elements = nelements(mesh)
    (; e_L, e_R) = cache
    if is_left
        # TODO: We cannot use `u[:, end, end]` here because for `PerElementFDSBP` `u` is a
        # `VectorOfArray` of vectors with different lengths, where `end` is not well-defined
        # and can give wrong results:
        # https://github.com/SciML/RecursiveArrayTools.jl/issues/454#issuecomment-2927845128
        return get_multiplied_node_vars(u, equations, e_R', :, N_elements)
    else
        return get_multiplied_node_vars(u, equations, e_L', :, 1)
    end
end

@inline function (::BoundaryConditionPeriodic)(u, x, t, mesh::OversetGridMesh, equations,
                                               solver, is_left, cache)
    u_left, u_right = u
    N_elements = nelements(mesh.mesh_right)
    _, solver_right = solver
    (; e_L, e_R) = cache
    if is_left
        # TODO: We cannot use `u_right[:, end, end]` here because for `PerElementFDSBP` `u` is a
        # `VectorOfArray` of vectors with different lengths, where `end` is not well-defined
        # and can give wrong results:
        # https://github.com/SciML/RecursiveArrayTools.jl/issues/454#issuecomment-2927845128
        return get_multiplied_node_vars(u_right, equations, e_R', :, N_elements)
    else
        return get_multiplied_node_vars(u_left, equations, e_L', :, 1)
    end
end

struct BoundaryConditionDoNothing end

"""
    boundary_condition_do_nothing = BoundaryConditionDoNothing()

Imposing no boundary condition just evaluates the flux at the inner state.
"""
const boundary_condition_do_nothing = BoundaryConditionDoNothing()

function Base.show(io::IO, ::BoundaryConditionDoNothing)
    print(io, "boundary_condition_do_nothing")
end

@inline function (::BoundaryConditionDoNothing)(u, x, t, mesh, equations, solver, is_left, cache)
    N_elements = nelements(mesh)
    (; e_L, e_R) = cache
    if is_left
        return get_multiplied_node_vars(u, equations, e_L', :, 1)
    else
        # TODO: We cannot use `u[:, end, end]` here because for `PerElementFDSBP` `u` is a
        # `VectorOfArray` of vectors with different lengths, where `end` is not well-defined
        # and can give wrong results:
        # https://github.com/SciML/RecursiveArrayTools.jl/issues/454#issuecomment-2927845128
        return get_multiplied_node_vars(u, equations, e_R', :, N_elements)
    end
end

@inline function (::BoundaryConditionDoNothing)(u, x, t, mesh::OversetGridMesh, equations,
                                                solver, is_left, cache)
    u_left, u_right = u
    N_elements = nelements(mesh.mesh_right)
    (; e_L, e_R) = cache
    if is_left
        return get_multiplied_node_vars(u_left, equations, e_L', :, 1)
    else
        # TODO: We cannot use `u_right[:, end, end]` here because for `PerElementFDSBP` `u` is a
        # `VectorOfArray` of vectors with different lengths, where `end` is not well-defined
        # and can give wrong results:
        # https://github.com/SciML/RecursiveArrayTools.jl/issues/454#issuecomment-2927845128
        return get_multiplied_node_vars(u_right, equations, e_R', :, N_elements)
    end
end

"""
    BoundaryConditionDirichlet(boundary_value_function)

Create a Dirichlet boundary condition that uses the function `boundary_value_function`
to specify the values at the boundary.
This can be used to create a boundary condition that specifies exact boundary values
by passing the exact solution of the equation.
The passed boundary value function will be called with the same arguments as an initial condition function is called, i.e., as
```julia
boundary_value_function(x, t, equations)
```
where `x` specifies the coordinates, `t` is the current time, and `equation` is the corresponding system of equations.

# Examples
```julia
julia> BoundaryConditionDirichlet(initial_condition_convergence_test)
```
"""
struct BoundaryConditionDirichlet{B}
    boundary_value_function::B
end

function Base.show(io::IO, ::BoundaryConditionDirichlet)
    print(io, "boundary_condition_dirichlet")
end
@inline function (boundary_condition::BoundaryConditionDirichlet)(u, x, t, mesh, equations,
                                                                  solver, is_left)
    return boundary_condition.boundary_value_function(x, t, equations)
end
