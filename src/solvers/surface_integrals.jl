abstract type AbstractSurfaceIntegral end

function calc_interface_flux!(surface_flux_values, u, mesh,
                              equations, integral::AbstractSurfaceIntegral, solver, cache)
    (; e_left, e_right) = cache
    for element in 2:(nelements(mesh) - 1)
        # left interface
        e_L = get_projection_operator(e_left, solver, element)
        e_R_m1 = get_projection_operator(e_right, solver, element - 1)
        u_ll = get_multiplied_node_vars(u, equations, e_R_m1', :, element - 1)
        u_rr = get_multiplied_node_vars(u, equations, e_L', :, element)

        f = integral.surface_flux(u_ll, u_rr, equations)
        set_node_vars!(surface_flux_values, f, equations, 1, element)
        set_node_vars!(surface_flux_values, f, equations, 2, element - 1)
        # right interface
        e_R = get_projection_operator(e_right, solver, element)
        e_L_p1 = get_projection_operator(e_left, solver, element + 1)
        u_ll = get_multiplied_node_vars(u, equations, e_R', :, element)
        u_rr = get_multiplied_node_vars(u, equations, e_L_p1', :, element + 1)
        f = integral.surface_flux(u_ll, u_rr, equations)
        set_node_vars!(surface_flux_values, f, equations, 2, element)
        set_node_vars!(surface_flux_values, f, equations, 1, element + 1)
    end
end

function calc_boundary_flux!(surface_flux_values, u, t, boundary_conditions, mesh,
                             equations, integral::AbstractSurfaceIntegral, solver, cache)
    (; x_neg, x_pos) = boundary_conditions
    (; e_left, e_right) = cache
    e_L = get_projection_operator(e_left, solver, 1)
    u_ll = x_neg(u, xmin(mesh), t, mesh, equations, solver, true, cache)
    u_rr = get_multiplied_node_vars(u, equations, e_L', :, 1)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 1, 1)

    e_R = get_projection_operator(e_right, solver, nelements(mesh))
    u_ll = get_multiplied_node_vars(u, equations, e_R', :, nelements(mesh))
    u_rr = x_pos(u, xmax(mesh), t, mesh, equations, solver, false, cache)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 2, nelements(mesh))
    return nothing
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
    return SurfaceIntegralStrongForm(surface_flux[1], surface_flux[2])
end
function SurfaceIntegralStrongForm(surface_flux)
    return SurfaceIntegralStrongForm(surface_flux, surface_flux)
end
SurfaceIntegralStrongForm() = SurfaceIntegralStrongForm(flux_central)

# This is M^{-1} * B * (f* - f) for `B = Diagonal([-1, 0, ..., 0, 1])` and `f* = [f_L^{num}, 0, ..., 0, f_R^{num}]`
# So basically a SAT.
function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::SurfaceIntegralStrongForm; left)
    M = mass_matrix(basis)
    boundary = left ? SummationByPartsOperators.xmin(basis) :
               SummationByPartsOperators.xmax(basis)
    boundary_projection = _interpolation_matrix([boundary], basis)[1, :]
    return M \ boundary_projection
end

function create_cache_surface_flux_values(mesh, equations, solver)
    surface_flux_values = zeros(real(solver), nvariables(equations), 2, nelements(mesh))
    return surface_flux_values
end

function create_cache(mesh, equations, solver, integral::SurfaceIntegralStrongForm)
    e_left, e_right = compute_projection_operators(solver)
    surface_operator_left = compute_integral_operator(mesh, solver, integral; left = true)
    surface_operator_right = compute_integral_operator(mesh, solver, integral; left = false)

    surface_flux_values = create_cache_surface_flux_values(mesh, equations, solver)
    return (; surface_operator_left, surface_operator_right, surface_flux_values, e_left,
            e_right)
end

@views function calc_surface_integral!(du, u, mesh, equations,
                                       ::SurfaceIntegralStrongForm, solver, cache)
    (; surface_operator_left, surface_operator_right, surface_flux_values, e_left, e_right) = cache
    for element in eachelement(mesh)
        e_L = get_projection_operator(e_left, solver, element)
        e_R = get_projection_operator(e_right, solver, element)
        u_L = get_multiplied_node_vars(u, equations, e_L', :, element)
        f_L = flux(u_L, equations)
        u_R = get_multiplied_node_vars(u, equations, e_R', :, element)
        f_R = flux(u_R, equations)
        surface_operator_left_ = get_integral_operator(surface_operator_left, solver,
                                                       element)
        surface_operator_right_ = get_integral_operator(surface_operator_right, solver,
                                                        element)
        for v in eachvariable(equations)
            du[v, :, element] .= du[v, :, element] +
                                 surface_operator_left_ *
                                 (surface_flux_values[v, 1, element] - f_L[v]) -
                                 surface_operator_right_ *
                                 (surface_flux_values[v, 2, element] - f_R[v])
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
    return SurfaceIntegralWeakForm(surface_flux[1], surface_flux[2])
end
SurfaceIntegralWeakForm(surface_flux) = SurfaceIntegralWeakForm(surface_flux, surface_flux)
SurfaceIntegralWeakForm() = SurfaceIntegralWeakForm(flux_central)

# This is M^{-1} * B * f* for `B = Diagonal([-1, 0, ..., 0, 1])` and `f* = [f_L^{num}, 0, ..., 0, f_R^{num}]`
function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::SurfaceIntegralWeakForm)
    M = mass_matrix(basis)
    x_l_ref = SummationByPartsOperators.xmin(basis)
    x_r_ref = SummationByPartsOperators.xmax(basis)
    R = _interpolation_matrix([x_l_ref, x_r_ref], basis)
    B = Diagonal([-1, 1])
    return -M \ (R' * B)
end

function create_cache(mesh, equations, solver, integral::SurfaceIntegralWeakForm)
    e_left, e_right = compute_projection_operators(solver)
    surface_operator = compute_integral_operator(mesh, solver, integral)

    surface_flux_values = create_cache_surface_flux_values(mesh, equations, solver)
    return (; surface_operator, surface_flux_values, e_left, e_right)
end

@views function calc_surface_integral!(du, u, mesh, equations,
                                       ::SurfaceIntegralWeakForm, solver, cache)
    (; surface_operator, surface_flux_values) = cache
    for element in eachelement(mesh)
        surface_operator_ = get_integral_operator(surface_operator, solver, element)
        for v in eachvariable(equations)
            du[v, :, element] .= du[v, :, element] +
                                 surface_operator_ * surface_flux_values[v, :, element]
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
