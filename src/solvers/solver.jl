include("volume_integrals.jl")
include("surface_integrals.jl")
include("dg.jl")
include("dgsem.jl")
include("per_element.jl")

# We define this function to allow for a different operator per element,
# This way we can use the same function for both classical `DG` and `PerElementFDSBP`,
# but without having to also store redundant information in the classical `DG` case.
get_basis(solver::DG, element) = solver.basis
get_basis(solver::PerElementFDSBP, element) = solver.basis.bases[element]

get_integral_operator(operator, solver, element) = operator
function get_integral_operator(operator, ::PerElementFDSBP, element)
    operator[element]
end

function compute_integral_operator(solver::DG, integral; kwargs...)
    compute_integral_operator(solver.basis, integral; kwargs...)
end
function compute_integral_operator(solver::PerElementFDSBP, integral; kwargs...)
    n_elements = length(solver.basis.bases)
    # Vector needed for `SurfaceIntegralStrongForm`, otherwise we need a `Matrix`.
    integral_operator = Vector{VecOrMat{real(solver)}}(undef, n_elements)
    for element in 1:n_elements
        integral_operator[element] = compute_integral_operator(get_basis(solver, element),
                                                               integral; kwargs...)
    end
    return integral_operator
end
