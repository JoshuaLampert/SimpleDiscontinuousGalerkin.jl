"""
    AbstractEquations{NDIMS, NVARS}

An abstract supertype of specific equations such as the linear advection equation.
The type parameters encode the number of spatial dimensions (`NDIMS`) and the
number of primary variables (`NVARS`) of the physics model.
"""
abstract type AbstractEquations{NDIMS, NVARS} end

@inline Base.ndims(::AbstractEquations{NDIMS}) where {NDIMS} = NDIMS

# Retrieve number of variables from equation instance
@inline nvariables(::AbstractEquations{NDIMS, NVARS}) where {NDIMS, NVARS} = NVARS

"""
    eachvariable(equations::AbstractEquations)

Return an iterator over the indices that specify the location in relevant data structures
for the variables in `equations`. In particular, not the variables themselves are returned.
"""
@inline eachvariable(equations::AbstractEquations) = Base.OneTo(nvariables(equations))

Base.broadcastable(equations::AbstractEquations) = (equations,)

"""
    get_name(equations::AbstractEquations)

Return the canonical, human-readable name for the given system of equations.
# Examples
```jldoctest
julia> SimpleDiscontinuousGalerkin.get_name(LinearAdvectionEquation1D(advection_velocity = 1.0))
"LinearAdvectionEquation1D"
```
"""
get_name(equations::AbstractEquations) = equations |> typeof |> nameof |> string

include("numerical_fluxes.jl")
include("linear_advection.jl")
