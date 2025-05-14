abstract type AbstractVolumeIntegral end

@views function calc_volume_integral!(du, u, mesh, equations,
                                      ::AbstractVolumeIntegral, dg, cache)
    (; volume_operator, f_all) = cache
    for element in eachelement(mesh)
        for node in eachnode(dg)
            f = flux(u[:, node, element], equations)
            for v in eachvariable(equations)
                f_all[v, node, element] = f[v]
            end
        end
        for v in eachvariable(equations)
            du[v, :, element] .= du[v, :, element] + volume_operator * f_all[v, :, element]
        end
    end
    return nothing
end

"""
    VolumeIntegralStrongForm()

The classical strong form volume integral type for FD/DG methods.
"""
struct VolumeIntegralStrongForm <: AbstractVolumeIntegral end

function create_cache(mesh, equations, solver, ::VolumeIntegralStrongForm)
    D = Matrix(solver.basis)
    volume_operator = -D
    f_all = zeros(real(solver), nvariables(equations), nnodes(solver), nelements(mesh))
    return (; volume_operator, f_all)
end

"""
    VolumeIntegralWeakForm()

The classical weak form volume integral type for DG methods as explained in
standard textbooks.

## References

- Kopriva (2009)
  Implementing Spectral Methods for Partial Differential Equations:
  Algorithms for Scientists and Engineers
  [doi: 10.1007/978-90-481-2261-5](https://doi.org/10.1007/978-90-481-2261-5)
- Hesthaven, Warburton (2007)
  Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and
  Applications
  [doi: 10.1007/978-0-387-72067-8](https://doi.org/10.1007/978-0-387-72067-8)

`VolumeIntegralWeakForm()` is only implemented for conserved terms as
non-conservative terms should always be discretized in conjunction with a flux-splitting scheme,
see [`VolumeIntegralFluxDifferencing`](@ref).
This treatment is required to achieve, e.g., entropy-stability or well-balancedness.
"""
struct VolumeIntegralWeakForm <: AbstractVolumeIntegral end

function create_cache(mesh, equations, solver, ::VolumeIntegralWeakForm)
    M = mass_matrix(solver.basis)
    D = Matrix(solver.basis)
    volume_operator = (M \ D') * M
    f_all = zeros(real(solver), nvariables(equations), nnodes(solver), nelements(mesh))
    return (; volume_operator, f_all)
end

abstract type AbstractSurfaceIntegral end

function calc_interface_flux!(surface_flux_values, u, mesh,
                              equations, integral::AbstractSurfaceIntegral, dg, cache)
    N = nnodes(dg)
    for element in 2:(nelements(mesh) - 1)
        # left interface
        u_ll = get_node_vars(u, equations, N, element - 1)
        u_rr = get_node_vars(u, equations, 1, element)
        f = integral.surface_flux(u_ll, u_rr, equations)
        set_node_vars!(surface_flux_values, f, equations, 1, element)
        set_node_vars!(surface_flux_values, f, equations, 2, element - 1)
        # right interface
        u_ll = get_node_vars(u, equations, N, element)
        u_rr = get_node_vars(u, equations, 1, element + 1)
        f = integral.surface_flux(u_ll, u_rr, equations)
        set_node_vars!(surface_flux_values, f, equations, 2, element)
        set_node_vars!(surface_flux_values, f, equations, 1, element + 1)
    end
end

function calc_boundary_flux!(surface_flux_values, u, t, boundary_conditions, mesh,
                             equations, integral::AbstractSurfaceIntegral, dg)
    (; x_neg, x_pos) = boundary_conditions
    u_ll = x_neg(u, xmin(mesh), t, equations, true)
    u_rr = get_node_vars(u, equations, 1, 1)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 1, 1)

    u_ll = get_node_vars(u, equations, nnodes(dg), nelements(mesh))
    u_rr = x_pos(u, xmax(mesh), t, equations, false)
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
struct SurfaceIntegralStrongForm{SurfaceFlux, SurfaceFluxBoundary} <: AbstractSurfaceIntegral
    surface_flux::SurfaceFlux
    surface_flux_boundary::SurfaceFluxBoundary
end

SurfaceIntegralStrongForm(surface_flux::Tuple) = SurfaceIntegralStrongForm(surface_flux[1], surface_flux[2])
SurfaceIntegralStrongForm(surface_flux) = SurfaceIntegralStrongForm(surface_flux, surface_flux)
SurfaceIntegralStrongForm() = SurfaceIntegralStrongForm(flux_central)

function create_cache(mesh, equations, solver, ::SurfaceIntegralStrongForm)
    M = mass_matrix(solver.basis)
    e_L = zeros(real(solver), nnodes(solver))
    e_L[1] = 1
    surface_operator_left = M \ e_L
    e_R = zeros(real(solver), nnodes(solver))
    e_R[end] = 1
    surface_operator_right = M \ e_R

    surface_flux_values = zeros(real(solver), nvariables(equations), 2, nelements(mesh))
    return (; surface_operator_left, surface_operator_right, surface_flux_values)
end

@views function calc_surface_integral!(du, u, mesh, equations,
                                       ::SurfaceIntegralStrongForm, dg, cache)
    (; surface_operator_left, surface_operator_right, surface_flux_values) = cache
    for element in eachelement(mesh)
        f_L = flux(u[:, 1, element], equations)
        f_R = flux(u[:, end, element], equations)
        for v in eachvariable(equations)
            du[v, :, element] .= du[v, :, element] +
                                 surface_operator_left *
                                 (surface_flux_values[v, 1, element] - f_L[v]) -
                                 surface_operator_right *
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
        print(io, "SurfaceIntegralStrongForm(", integral.surface_flux, ")")
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

SurfaceIntegralWeakForm(surface_flux::Tuple) = SurfaceIntegralWeakForm(surface_flux[1], surface_flux[2])
SurfaceIntegralWeakForm(surface_flux) = SurfaceIntegralWeakForm(surface_flux, surface_flux)
SurfaceIntegralWeakForm() = SurfaceIntegralWeakForm(flux_central)

function create_cache(mesh, equations, solver, ::SurfaceIntegralWeakForm)
    M = mass_matrix(solver.basis)
    R = zeros(real(solver), 2, nnodes(solver))
    R[1, 1] = R[end, end] = 1
    B = Diagonal([-1, 1])
    surface_operator = -M \ (R' * B)

    surface_flux_values = zeros(real(solver), nvariables(equations), 2, nelements(mesh))
    return (; surface_operator, surface_flux_values)
end

@views function calc_surface_integral!(du, u, mesh, equations,
                                       ::SurfaceIntegralWeakForm, dg, cache)
    (; surface_operator, surface_flux_values) = cache
    for v in eachvariable(equations)
        du[v, :, :] .= du[v, :, :] + surface_operator * surface_flux_values[v, :, :]
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", integral::SurfaceIntegralWeakForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "SurfaceIntegralWeakForm(", integral.surface_flux, ")")
    end
end
