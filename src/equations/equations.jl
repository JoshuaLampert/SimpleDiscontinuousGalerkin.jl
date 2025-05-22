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

"""
    cons2cons(u, equations)

Return the conservative variables `u`. While this function is as trivial as `identity`,
it is also as useful.
"""
@inline cons2cons(u, ::AbstractEquations) = u

# Add methods to show some information on systems of equations.
function Base.show(io::IO, equations::AbstractEquations)
    # Since this is not performance-critical, we can use `@nospecialize` to reduce latency.
    @nospecialize equations # reduce precompilation time

    print(io, get_name(equations), " with ")
    if nvariables(equations) == 1
        print(io, "one variable")
    else
        print(io, nvariables(equations), " variables")
    end
end

function Base.show(io::IO, ::MIME"text/plain", equations::AbstractEquations)
    # Since this is not performance-critical, we can use `@nospecialize` to reduce latency.
    @nospecialize equations # reduce precompilation time

    if get(io, :compact, false)
        show(io, equations)
    else
        println(io, get_name(equations))
        print(io, "#variables: ", nvariables(equations))
        for variable in eachvariable(equations)
            println()
            print("    variable " * string(variable), ": ",
                  varnames(cons2cons, equations)[variable])
        end
    end
end

include("numerical_fluxes.jl")
include("linear_advection.jl")
