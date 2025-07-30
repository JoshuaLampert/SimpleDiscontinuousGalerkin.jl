abstract type AbstractSurfaceIntegral end

function calc_interface_flux!(surface_flux_values, u, mesh,
                              equations, integral::AbstractSurfaceIntegral, solver, cache)
    for element in 2:(nelements(mesh) - 1)
        N_element = nnodes(solver, element)
        N_element_m1 = nnodes(solver, element - 1)
        # left interface
        u_ll = get_node_vars(u, equations, N_element_m1, element - 1)
        u_rr = get_node_vars(u, equations, 1, element)
        f = integral.surface_flux(u_ll, u_rr, equations)
        set_node_vars!(surface_flux_values, f, equations, 1, element)
        set_node_vars!(surface_flux_values, f, equations, 2, element - 1)
        # right interface
        u_ll = get_node_vars(u, equations, N_element, element)
        u_rr = get_node_vars(u, equations, 1, element + 1)
        f = integral.surface_flux(u_ll, u_rr, equations)
        set_node_vars!(surface_flux_values, f, equations, 2, element)
        set_node_vars!(surface_flux_values, f, equations, 1, element + 1)
    end
end

function calc_boundary_flux!(surface_flux_values, u, t, boundary_conditions, mesh,
                             equations, integral::AbstractSurfaceIntegral, solver)
    (; x_neg, x_pos) = boundary_conditions
    u_ll = x_neg(u, xmin(mesh), t, mesh, equations, solver, true)
    u_rr = get_node_vars(u, equations, 1, 1)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 1, 1)

    u_ll = get_node_vars(u, equations, nnodes(solver, nelements(mesh)), nelements(mesh))
    u_rr = x_pos(u, xmax(mesh), t, mesh, equations, solver, false)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 2, nelements(mesh))
    return nothing
end

function calc_surface_integral!(du, u, mesh, equations, integral::SurfaceIntegralWeakForm,
                                solver, cache)
    calc_surface_integral!(du, u, mesh, equations, integral, solver, cache.surface_operator,
                           cache.surface_flux_values)
end

function calc_surface_integral!(du, u, mesh, equations, integral::SurfaceIntegralStrongForm,
                                solver, cache)
    calc_surface_integral!(du, u, mesh, equations, integral, solver,
                           cache.surface_operator_left, cache.surface_operator_right,
                           cache.surface_flux_values)
end

"""
    SurfaceIntegralStrongForm(surface_flux=flux_central, surface_flux_boundary=surface_flux)

The classical strong form surface integral type for FD/DG methods.
It uses `surface_flux` for the interior fluxes and `surface_flux_boundary` for the
boundary fluxes.

See also [`VolumeIntegralStrongForm`](@ref).
"""
struct SurfaceIntegralStrongForm{SurfaceFlux, SurfaceFluxBoundary} <:
       AbstractSurfaceIntegral
    surface_flux::SurfaceFlux
    surface_flux_boundary::SurfaceFluxBoundary
end

function SurfaceIntegralStrongForm(surface_flux::Tuple)
    SurfaceIntegralStrongForm(surface_flux[1], surface_flux[2])
end
function SurfaceIntegralStrongForm(surface_flux)
    SurfaceIntegralStrongForm(surface_flux, surface_flux)
end
SurfaceIntegralStrongForm() = SurfaceIntegralStrongForm(flux_central)

# This is M^{-1} * B * (f* - f) for `B = Diagonal([-1, 0, ..., 0, 1])` and `f* = [f_L^{num}, 0, ..., 0, f_R^{num}]`
# So basically a SAT.
# This assumes Lobatto-type nodes, where the first and last node are the boundaries.
function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::SurfaceIntegralStrongForm; left)
    M = mass_matrix(basis)
    unit_vector = zeros(eltype(M), length(grid(basis)))
    left ? unit_vector[begin] = 1 : unit_vector[end] = 1
    return M \ unit_vector
end

function create_cache_surface_flux_values(mesh, equations, solver)
    surface_flux_values = zeros(real(solver), nvariables(equations), 2, nelements(mesh))
    return surface_flux_values
end

function create_cache(mesh, equations, solver, integral::SurfaceIntegralStrongForm)
    surface_operator_left = compute_integral_operator(mesh, solver, integral; left = true)
    surface_operator_right = compute_integral_operator(mesh, solver, integral; left = false)

    surface_flux_values = create_cache_surface_flux_values(mesh, equations, solver)
    return (; surface_operator_left, surface_operator_right, surface_flux_values)
end

# TODO: Here, we would like to use `@views` to avoid allocations, but there is currently
# a bug in RecursiveArrayTools.jl: https://github.com/SciML/RecursiveArrayTools.jl/issues/453
function calc_surface_integral!(du, u, mesh, equations,
                                ::SurfaceIntegralStrongForm, solver, surface_operator_left,
                                surface_operator_right, surface_flux_values)
    for element in eachelement(mesh)
        f_L = flux(u[:, 1, element], equations)
        # TODO: We cannot use `u[:, end, element]` here because for `PerElementFDSBP` `u` is a
        # `VectorOfArray` of vectors with different lengths, where `end` is not well-defined
        # and can give wrong results:
        # https://github.com/SciML/RecursiveArrayTools.jl/issues/454#issuecomment-2927845128
        N = nnodes(solver, element)
        f_R = flux(u[:, N, element], equations)
        surface_operator_left_ = get_integral_operator(surface_operator_left, solver,
                                                       element)
        surface_operator_right_ = get_integral_operator(surface_operator_right, solver,
                                                        element)
        for v in eachvariable(equations)
            # TODO: We would like to use broadcasting here:
            # du[v, :, element] .= du[v, :, element] +
            #                      surface_operator_left_ *
            #                     (surface_flux_values[v, 1, element] - f_L[v]) -
            #                     surface_operator_right_ *
            #                     (surface_flux_values[v, 2, element] - f_R[v])
            # but there are currently issues with RecursiveArrayTools.jl:
            # https://github.com/SciML/RecursiveArrayTools.jl/issues/453 and https://github.com/SciML/RecursiveArrayTools.jl/issues/454
            du_update = surface_operator_left_ *
                        (surface_flux_values[v, 1, element] - f_L[v]) -
                        surface_operator_right_ *
                        (surface_flux_values[v, 2, element] - f_R[v])
            for node in eachnode(solver, element)
                du[v, node, element] += du_update[node]
            end
        end
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", integral::SurfaceIntegralStrongForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "SurfaceIntegralStrongForm(", integral.surface_flux, ", ",
              integral.surface_flux_boundary, ")")
    end
end

"""
    SurfaceIntegralWeakForm(surface_flux=flux_central, surface_flux_boundary=surface_flux)

The classical weak form surface integral type for DG methods as explained in standard
textbooks. It uses `surface_flux` for the interior fluxes and `surface_flux_boundary` for the
boundary fluxes.

See also [`VolumeIntegralWeakForm`](@ref).

## References

- Kopriva (2009)
  Implementing Spectral Methods for Partial Differential Equations:
  Algorithms for Scientists and Engineers
  [doi: 10.1007/978-90-481-2261-5](https://doi.org/10.1007/978-90-481-2261-5)
- Hesthaven, Warburton (2007)
  Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and
  Applications
  [doi: 10.1007/978-0-387-72067-8](https://doi.org/10.1007/978-0-387-72067-8)
"""
struct SurfaceIntegralWeakForm{SurfaceFlux, SurfaceFluxBoundary} <: AbstractSurfaceIntegral
    surface_flux::SurfaceFlux
    surface_flux_boundary::SurfaceFluxBoundary
end

function SurfaceIntegralWeakForm(surface_flux::Tuple)
    SurfaceIntegralWeakForm(surface_flux[1], surface_flux[2])
end
SurfaceIntegralWeakForm(surface_flux) = SurfaceIntegralWeakForm(surface_flux, surface_flux)
SurfaceIntegralWeakForm() = SurfaceIntegralWeakForm(flux_central)

# This is M^{-1} * B * f* for `B = Diagonal([-1, 0, ..., 0, 1])` and `f* = [f_L^{num}, 0, ..., 0, f_R^{num}]`
# This assumes Lobatto-type nodes, where the first and last node are the boundaries.
function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::SurfaceIntegralWeakForm)
    M = mass_matrix(basis)
    R = zeros(eltype(M), 2, length(grid(basis)))
    R[1, 1] = R[end, end] = 1
    B = Diagonal([-1, 1])
    return -M \ (R' * B)
end

function create_cache(mesh, equations, solver, integral::SurfaceIntegralWeakForm)
    surface_operator = compute_integral_operator(mesh, solver, integral)

    surface_flux_values = create_cache_surface_flux_values(mesh, equations, solver)
    return (; surface_operator, surface_flux_values)
end

# TODO: Here, we would like to use `@views` to avoid allocations, but there is currently
# a bug in RecursiveArrayTools.jl: https://github.com/SciML/RecursiveArrayTools.jl/issues/453
function calc_surface_integral!(du, u, mesh, equations,
                                ::SurfaceIntegralWeakForm, solver, surface_operator,
                                surface_flux_values)
    for element in eachelement(mesh)
        surface_operator_ = get_integral_operator(surface_operator, solver, element)
        for v in eachvariable(equations)
            # TODO: We would like to use broadcasting here:
            # du[v, :, element] .= du[v, :, element] +
            #                      surface_operator_ * surface_flux_values[v, :, element]
            # but there are currently issues with RecursiveArrayTools.jl:
            # https://github.com/SciML/RecursiveArrayTools.jl/issues/453 and https://github.com/SciML/RecursiveArrayTools.jl/issues/454
            du_update = surface_operator_ * surface_flux_values[v, :, element]
            for node in eachnode(solver, element)
                du[v, node, element] += du_update[node]
            end
        end
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", integral::SurfaceIntegralWeakForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "SurfaceIntegralWeakForm(", integral.surface_flux, ", ",
              integral.surface_flux_boundary, ")")
    end
end
