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

function compute_integral_operator(mesh, solver::DG, integral; kwargs...)
    compute_integral_operator(solver.basis, integral; kwargs...)
end
function compute_integral_operator(mesh, solver::PerElementFDSBP, integral; kwargs...)
    n_elements = nelements(mesh)
    # Vector needed for `SurfaceIntegralStrongForm`, otherwise we need a `Matrix`.
    integral_operator = Vector{VecOrMat{real(solver)}}(undef, n_elements)
    for element in eachelement(mesh)
        integral_operator[element] = compute_integral_operator(get_basis(solver, element),
                                                               integral; kwargs...)
    end
    return integral_operator
end

function calc_error_norms(u, t, initial_condition, mesh, equations,
                          solver, cache)
    u_exact = similar(u)
    compute_coefficients!(u_exact, initial_condition, t, mesh, equations, solver, cache)
    l2_error = zeros(real(solver), nvariables(equations))
    linf_error = zeros(real(solver), nvariables(equations))
    for v in eachvariable(equations)
        l2_error[v] = zero(real(solver))
        linf_error[v] = zero(real(solver))
        for element in eachelement(mesh)
            # TODO: Broadcast (issue in RecursiveArrayTools.jl)
            diff = zeros(real(solver), nnodes(solver, element))
            for i in eachnode(solver, element)
                diff[i] = u[v, i, element] - u_exact[v, i, element]
            end
            l2_error[v] += cache.jacobian[element] *
                           integrate(abs2, diff, get_basis(solver, element))
            linf_error[v] = max(linf_error[v], maximum(abs.(diff)))
        end
        l2_error[v] = sqrt(l2_error[v])
    end
    return l2_error, linf_error
end

function calc_sources!(du, u, t, source_terms::Nothing, mesh, equations, solver, cache)
    return nothing
end

function calc_sources!(du, u, t, source_terms, mesh, equations, solver, cache)
    node_coordinates = cache.node_coordinates

    for element in eachelement(mesh)
        for i in eachnode(solver, element)
            u_local = get_node_vars(u, equations, i, element)
            x_local = get_node_coords(node_coordinates, equations, solver,
                                      i, element)
            du_local = source_terms(u_local, x_local, t, equations)
            for v in eachvariable(equations)
                du[v, i, element] = du[v, i, element] + du_local[v]
            end
        end
    end

    return nothing
end
