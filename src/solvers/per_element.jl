"""
    PerElementBasis{BasisType}

Basis, which can hold a different SBP operator for each element.
This is used in [`PerElementFDSBP`](@ref) to allow for different bases on each element.
"""
struct PerElementBasis{BasisType}
    bases::Vector{BasisType}
end

"""
    PerElementFDSBP(bases::Vector{BasisType};
                    surface_flux = flux_central,
                    surface_integral = SurfaceIntegralWeakForm(surface_flux),
                    volume_integral = VolumeIntegralWeakForm()) where BasisType

Create a discontinuous Galerkin method using different bases for each element.
This is like [`FDSBP`](@ref), but allows for a different SBP operator on each element.
See also: [`PerElementBasis`](@ref).
"""
const PerElementFDSBP = DG{Basis} where {Basis <: PerElementBasis}
function PerElementFDSBP(bases::Vector{BasisType};
                         surface_flux = flux_central,
                         surface_integral = SurfaceIntegralWeakForm(surface_flux),
                         volume_integral = VolumeIntegralWeakForm()) where {BasisType}
    return DG{PerElementBasis{typeof(bases)}, typeof(surface_integral),
              typeof(volume_integral)}(PerElementBasis(bases), surface_integral,
                                       volume_integral)
end

function Base.summary(io::IO, dg::PerElementFDSBP)
    print(io, "PerElementFDSBP(bases=$(dg.basis.bases))")
end

function create_cache(mesh, equations, dg::PerElementFDSBP, initial_condition,
                      boundary_conditions)
    @assert length(dg.basis.bases)==nelements(mesh) "Number of bases must match number of elements in the mesh"
    dx = element_spacing(mesh) # length of each element
    # We need a `Vector{Vector}` to account for potentially different number of nodes for each element
    # compute all mapped nodes
    x = Vector{Vector{real(dg)}}(undef, nelements(mesh))
    jacobian = zeros(real(dg), nelements(mesh))
    for element in eachelement(mesh)
        x_l = xmin(mesh) + (element - 1) * dx
        D = get_basis(dg, element)
        nodes_basis = grid(D)
        dx_basis = last(nodes_basis) - first(nodes_basis)
        jacobian[element] = dx_u / dx_basis
        x[element] = zeros(length(nodes_basis))
        for j in eachindex(nodes_basis)
            x[j, element] = x_l + jacobian * (nodes_basis[j] - first(nodes_basis))
        end
    end
    cache = (; jacobian, node_coordinates = x,
             create_cache(mesh, equations, dg, dg.volume_integral)...,
             create_cache(mesh, equations, dg, dg.surface_integral)...)
    return cache
end

function apply_jacobian!(du, mesh, equations, dg::PerElementFDSBP, cache)
    (; jacobian) = cache
    for element in eachelement(mesh)
        @. du[element] = du[element] / jacobian[element]
    end
end
