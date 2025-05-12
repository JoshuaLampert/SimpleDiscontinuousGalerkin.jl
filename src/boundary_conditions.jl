function digest_boundary_conditions(boundary_conditions, mesh::Mesh{1}, solver, cache)
    (; x_neg = boundary_conditions, x_pos = boundary_conditions)
end

function digest_boundary_conditions(boundary_conditions::NTuple{2, Any}, mesh::Mesh{1},
                                    solver, cache)
    (; x_neg = boundary_conditions[1], x_pos = boundary_conditions[2])
end
function digest_boundary_conditions(boundary_conditions::NamedTuple{Keys, ValueTypes},
                                    mesh::Mesh{1}, solver,
                                    cache) where {Keys, ValueTypes <: NTuple{2, Any}}
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

@inline function (::BoundaryConditionPeriodic)(u, t, equations, is_left)
    if is_left
        return u[:, end, end]
    else
        return u[:, 1, 1]
    end
end
