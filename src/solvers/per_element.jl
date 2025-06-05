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

function Base.summary(io::IO, solver::PerElementFDSBP)
    print(io, "PerElementFDSBP(bases=$(solver.basis.bases))")
end

grid(solver::PerElementFDSBP, element) = grid(solver.basis, element)

function create_cache(mesh, equations, solver::PerElementFDSBP, initial_condition,
                      boundary_conditions)
    @assert length(solver.basis.bases)==nelements(mesh) "Number of bases must match number of elements in the mesh"
    # We need a `Vector{Vector}` to account for potentially different number of nodes for each element
    # compute all mapped nodes
    x = VectorOfArray([zeros(real(solver), nnodes(solver, element))
                       for element in eachelement(mesh)])
    jacobian = zeros(real(solver), nelements(mesh))
    for element in eachelement(mesh)
        x_l = left_element_boundary(mesh, element)
        D = get_basis(solver, element)
        nodes_basis = grid(D)
        dx_basis = last(nodes_basis) - first(nodes_basis)
        dx = element_spacing(mesh, element) # length of the element
        jacobian[element] = dx / dx_basis
        for j in eachindex(nodes_basis)
            x[j, element] = x_l + jacobian[element] * (nodes_basis[j] - first(nodes_basis))
        end
    end
    cache = (; jacobian, node_coordinates = x,
             create_cache(mesh, equations, solver, solver.volume_integral)...,
             create_cache(mesh, equations, solver, solver.surface_integral)...)
    return cache
end

# We need different data structure because we have a different number of nodes per element
# `get_node_coords`, `get_node_vars`, and `set_node_vars!` can be reused because we can
# access the `Array{RealT, 3}` in the same way as the `VectorOfArray` structure
# (`v` first, node variable(s) in the middle, `element` last).
function allocate_coefficients(mesh::AbstractMesh, equations, solver::PerElementFDSBP)
    u = [zeros(real(solver), nvariables(equations), nnodes(solver, element))
         for element in eachelement(mesh)]
    return VectorOfArray(u)
end

function get_variable(u, v, ::PerElementFDSBP)
    return collect(Iterators.flatten(getindex.(parent(u), v, :)))
end
