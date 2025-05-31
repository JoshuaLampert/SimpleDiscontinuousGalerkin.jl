"""
    PerElementBasis{BasisType}

Basis, which can hold a different SBP operator for each element.
This is used in [`PerElementFDSBP`](@ref) to allow for different bases on each element.
"""
struct PerElementBasis{BasisType}
    bases::Vector{BasisType}
end

@inline Base.real(basis::PerElementBasis) = real(first(basis.bases))
grid(basis::PerElementBasis, element) = grid(basis.bases[element])

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
    return DG{PerElementBasis{eltype(bases)}, typeof(surface_integral),
              typeof(volume_integral)}(PerElementBasis(bases), surface_integral,
                                       volume_integral)
end

function Base.summary(io::IO, dg::PerElementFDSBP)
    print(io, "PerElementFDSBP(bases=$(dg.basis.bases))")
end

grid(solver::PerElementFDSBP, element) = grid(solver.basis, element)

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
        jacobian[element] = dx / dx_basis
        x[element] = zeros(real(dg), length(nodes_basis))
        for j in eachindex(nodes_basis)
            x[element][j] = x_l + jacobian[element] * (nodes_basis[j] - first(nodes_basis))
        end
    end
    cache = (; jacobian, node_coordinates = VectorOfArray(x),
             create_cache(mesh, equations, dg, dg.volume_integral)...,
             create_cache(mesh, equations, dg, dg.surface_integral)...)
    return cache
end

function apply_jacobian!(du, mesh, equations, dg::PerElementFDSBP, cache)
    (; jacobian) = cache
    for element in eachelement(mesh)
        for i in eachnode(dg, element)
            for v in eachvariable(equations)
                du[v, i, element] = du[v, i, element] / jacobian[element]
            end
        end
    end
end

# We need different data structure because we have a different number of nodes per element
# `get_node_coords`, `get_node_vars`, and `set_node_vars!` can be reused because we can
# access the `Array{RealT, 3}` in the same way as the `VectorOfArray` structure
# (`v` first, node variable(s) in the middle, `element` last).
function allocate_coefficients(mesh::Mesh, equations, solver::PerElementFDSBP)
    u = [zeros(real(solver), nvariables(equations), nnodes(solver, element))
         for element in eachelement(mesh)]
    return VectorOfArray(u)
end

function get_variable(u, v, ::PerElementFDSBP)
    # TODO: Get only `v`
    return collect(Iterators.flatten(parent(u)))
end
