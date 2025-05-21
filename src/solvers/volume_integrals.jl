abstract type AbstractVolumeIntegral end

@views function calc_volume_integral!(du, u, mesh, equations,
                                      ::AbstractVolumeIntegral, dg, cache)
    (; volume_operator, f_all) = cache
    for element in eachelement(mesh)
        for node in eachnode(dg)
            u_node = get_node_vars(u, equations, node, element)
            f = flux(u_node, equations)
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

"""
    VolumeIntegralFluxDifferencing(volume_flux=flux_central)

Volume integral type for DG methods based on SBP operators and flux differencing
using a symmetric two-point `volume_flux`. This `volume_flux` needs to satisfy
the interface of numerical fluxes.

To be used together with [`SurfaceIntegralWeakForm`](@ref).
"""
struct VolumeIntegralFluxDifferencing{VolumeFlux} <:
       AbstractVolumeIntegral
    volume_flux::VolumeFlux
end

VolumeIntegralFluxDifferencing() = VolumeIntegralFluxDifferencing(flux_central)

function create_cache(mesh, equations, solver, ::VolumeIntegralFluxDifferencing)
    weights = diag(mass_matrix(solver.basis))
    D = Matrix(solver.basis)
    D_split = 2 * D
    D_split[1, 1] += 1 / weights[1]
    D_split[end, end] -= 1 / weights[end]
    return (; D_split)
end

# Subtract D_split * f*
@views function calc_volume_integral!(du, u, mesh, equations,
                                      integral::VolumeIntegralFluxDifferencing, dg, cache)
    (; D_split) = cache
    for element in eachelement(mesh)
        for node in eachnode(dg)
            u_node = get_node_vars(u, equations, node, element)
            for node_2 in eachnode(dg)
                u_node_2 = get_node_vars(u, equations, node_2, element)
                f = integral.volume_flux(u_node, u_node_2, equations)
                for v in eachvariable(equations)
                    du[v, node, element] = du[v, node, element] -
                                           D_split[node, node_2] * f[v]
                end
            end
        end
    end
    return nothing
end
