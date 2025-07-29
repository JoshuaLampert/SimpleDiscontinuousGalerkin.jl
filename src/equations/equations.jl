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
julia> SimpleDiscontinuousGalerkin.get_name(LinearAdvectionEquation1D(1.0))
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

"""
    default_analysis_errors(equations)

Default analysis errors used by the [`AnalysisCallback`](@ref).
"""
default_analysis_errors(::AbstractEquations) = (:l2_error, :linf_error)

"""
    default_analysis_integrals(equations)

Default analysis integrals used by the [`AnalysisCallback`](@ref).
"""
default_analysis_integrals(::AbstractEquations) = (entropy, entropy_timederivative)

function default_analysis_integrals(::AbstractEquations{1, 1})
    (mass, entropy, entropy_timederivative)
end

"""
    mass(u, equations)

Return the mass of the conservative variables `u` for the given system of equations.
"""
function mass(u, ::AbstractEquations{1, 1})
    return u[1]
end

"""
    entropy(u, equations)

Return the entropy of the conservative variables `u` for the given system of equations.
"""
function entropy(u, ::AbstractEquations{1, 1})
    return 0.5 * u[1]^2
end

"""
    entropy_timederivative

The semi-discrete time derivative of the entropy, which can be used in the
[`AnalysisCallback`](@ref) to compute the time derivative of the entropy.
"""
function entropy_timederivative end

"""
    cons2entropy(u, equations)

Return the entropy variables from the conservative variables `u` for the given system of equations.
The entropy variables are defined as the derivative of the entropy with respect to the conservative variables.
"""
cons2entropy(u, equations::AbstractEquations{1}) = u

include("numerical_fluxes.jl")
include("linear_advection.jl")
include("Maxwell.jl")
