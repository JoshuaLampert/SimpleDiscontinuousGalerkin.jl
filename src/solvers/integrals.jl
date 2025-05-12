abstract type AbstractVolumeIntegral end

"""
    VolumeIntegralStrongForm()

The classical strong form volume integral type for FD/DG methods.
"""
struct VolumeIntegralStrongForm <: AbstractVolumeIntegral end

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
    return (; volume_operator)
end

function calc_volume_integral!(du, u, mesh, equations,
                               ::VolumeIntegralWeakForm, dg, cache)
    (; volume_operator) = cache
    for v in eachvariable(equations)
        du[v, :, :] .= du[v, :, :] + volume_operator * u[v, :, :]
    end
    return nothing
end

abstract type AbstractSurfaceIntegral end

"""
    SurfaceIntegralWeakForm(surface_flux=flux_central)

The classical weak form surface integral type for DG methods as explained in standard
textbooks.

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
struct SurfaceIntegralWeakForm{SurfaceFlux} <: AbstractSurfaceIntegral
    surface_flux::SurfaceFlux
end

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

function calc_interface_flux!(surface_flux_values, u, mesh,
                              equations, integral::SurfaceIntegralWeakForm, dg, cache)
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
                             equations, integral::SurfaceIntegralWeakForm, dg)
    (; x_neg, x_pos) = boundary_conditions
    u_ll = x_neg(u, t, equations, true)
    u_rr = get_node_vars(u, equations, 1, 1)
    f = integral.surface_flux(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 1, 1)

    u_ll = get_node_vars(u, equations, nnodes(dg), nelements(mesh))
    u_rr = x_pos(u, t, equations, false)
    f = integral.surface_flux(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values, f, equations, 2, nelements(mesh))
    return nothing
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

"""
    SurfaceIntegralStrongForm(surface_flux=flux_central)

The classical strong form surface integral type for FD/DG methods.

See also [`VolumeIntegralStrongForm`](@ref).
"""
struct SurfaceIntegralStrongForm{SurfaceFlux} <: AbstractSurfaceIntegral
    surface_flux::SurfaceFlux
end

SurfaceIntegralStrongForm() = SurfaceIntegralStrongForm(flux_central)

function Base.show(io::IO, ::MIME"text/plain", integral::SurfaceIntegralStrongForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "SurfaceIntegralStrongForm(", integral.surface_flux, ")")
    end
end
