abstract type AbstractVolumeIntegral end

"""
    VolumeIntegralStrongForm()

The classical strong form volume integral type for FD/DG methods.
"""
struct VolumeIntegralStrongForm <: AbstractVolumeIntegral end

function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::VolumeIntegralStrongForm)
    D = Matrix(basis)
    return -D
end

function Base.show(io::IO, ::MIME"text/plain", integral::VolumeIntegralStrongForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "VolumeIntegralStrongForm")
    end
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
"""
struct VolumeIntegralWeakForm <: AbstractVolumeIntegral end

function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::VolumeIntegralWeakForm)
    M = mass_matrix(basis)
    D = Matrix(basis)
    return (M \ D') * M
end

function Base.show(io::IO, ::MIME"text/plain", integral::VolumeIntegralWeakForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "VolumeIntegralWeakForm")
    end
end

function create_cache(mesh, equations, solver,
                      integral::Union{VolumeIntegralStrongForm,
                                      VolumeIntegralWeakForm})
    volume_operator = compute_integral_operator(solver, integral)
    f_all = allocate_coefficients(mesh, equations, solver)
    return (; volume_operator, f_all)
end

# TODO: Here, we would like to use `@views` to avoid allocations, but there is currently
# a bug in RecursiveArrayTools.jl: https://github.com/SciML/RecursiveArrayTools.jl/issues/453
function calc_volume_integral!(du, u, mesh, equations,
                                      ::Union{VolumeIntegralStrongForm,
                                              VolumeIntegralWeakForm},
                                      dg, cache)
    (; volume_operator, f_all) = cache
    for element in eachelement(mesh)
        volume_operator_ = get_integral_operator(volume_operator, dg, element)
        for node in eachnode(dg, element)
            u_node = get_node_vars(u, equations, node, element)
            f = flux(u_node, equations)
            for v in eachvariable(equations)
                f_all[v, node, element] = f[v]
            end
        end
        for v in eachvariable(equations)
            # TODO: We would like to use broadcasting here:
            # du[v, :, element] .= du[v, :, element] + volume_operator_ * f_all[v, :, element]
            # but there are currently issues with RecursiveArrayTools.jl:
            # https://github.com/SciML/RecursiveArrayTools.jl/issues/453 and https://github.com/SciML/RecursiveArrayTools.jl/issues/454
            du_update = volume_operator_ * f_all[v, :, element]
            for node in eachnode(dg, element)
                du[v, node, element] += du_update[node]
            end
        end
    end
    return nothing
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

function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::VolumeIntegralFluxDifferencing)
    weights = diag(mass_matrix(basis))
    D = Matrix(basis)
    D_split = 2 * D
    D_split[1, 1] += 1 / weights[1]
    D_split[end, end] -= 1 / weights[end]
    return D_split
end

function Base.show(io::IO, ::MIME"text/plain", integral::VolumeIntegralFluxDifferencing)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "VolumeIntegralFluxDifferencing(", integral.volume_flux, ")")
    end
end

"""
    VolumeIntegralFluxDifferencingStrongForm(volume_flux=flux_central)

Volume integral type for DG methods based on SBP operators and flux differencing
using a symmetric two-point `volume_flux`. This `volume_flux` needs to satisfy
the interface of numerical fluxes.

This is the strong formulation, which means it should be used together with [`SurfaceIntegralStrongForm`](@ref).
"""
struct VolumeIntegralFluxDifferencingStrongForm{VolumeFlux} <:
       AbstractVolumeIntegral
    volume_flux::VolumeFlux
end

function VolumeIntegralFluxDifferencingStrongForm()
    VolumeIntegralFluxDifferencingStrongForm(flux_central)
end

function compute_integral_operator(basis::AbstractDerivativeOperator,
                                   ::VolumeIntegralFluxDifferencingStrongForm)
    D = Matrix(basis)
    return 2 * D
end

function Base.show(io::IO, ::MIME"text/plain",
                   integral::VolumeIntegralFluxDifferencingStrongForm)
    @nospecialize integral # reduce precompilation time

    if get(io, :compact, false)
        show(io, integral)
    else
        print(io, "VolumeIntegralFluxDifferencingStrongForm(", integral.volume_flux, ")")
    end
end

function create_cache(mesh, equations, solver,
                      integral::Union{VolumeIntegralFluxDifferencing,
                                      VolumeIntegralFluxDifferencingStrongForm})
    volume_operator = compute_integral_operator(solver, integral)
    return (; volume_operator)
end

# Subtract D_split * f^{vol} for `VolumeIntegralFluxDifferencing` and 2 * D * f^{vol} for
# `VolumeIntegralFluxDifferencingStrongForm`.
function calc_volume_integral!(du, u, mesh, equations,
                                      integral::Union{VolumeIntegralFluxDifferencing,
                                                      VolumeIntegralFluxDifferencingStrongForm},
                                      dg, cache)
    (; volume_operator) = cache
    for element in eachelement(mesh)
        volume_operator_ = get_integral_operator(volume_operator, dg, element)
        for node in eachnode(dg, element)
            u_node = get_node_vars(u, equations, node, element)
            for node_2 in eachnode(dg, element)
                u_node_2 = get_node_vars(u, equations, node_2, element)
                f_vol = integral.volume_flux(u_node, u_node_2, equations)
                for v in eachvariable(equations)
                    du[v, node, element] = du[v, node, element] -
                                           volume_operator_[node, node_2] * f_vol[v]
                end
            end
        end
    end
    return nothing
end
