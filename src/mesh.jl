abstract type AbstractMesh{NDIMS, RealT} end

@inline Base.ndims(::AbstractMesh{NDIMS}) where {NDIMS} = NDIMS
@inline eachelement(mesh::AbstractMesh) = Base.OneTo(nelements(mesh))
@inline Base.real(::AbstractMesh{NDIMS, RealT}) where {NDIMS, RealT} = RealT

function Base.show(io::IO, mesh::AbstractMesh{RealT}) where {RealT}
    print(io, nameof(typeof(mesh)), "{", RealT, "} with ", nelements(mesh), " elements")
    print(io, " from ", xmin(mesh), " to ", xmax(mesh))
end

function Base.show(io::IO, ::MIME"text/plain",
                   mesh::AbstractMesh{NDIMS, RealT}) where {NDIMS, RealT}
    if get(io, :compact, false)
        show(io, mesh)
    else
        println(io, string(nameof(typeof(mesh))), "{", RealT, "} ")
        println(io, "    xmin: ", xmin(mesh))
        println(io, "    xmax: ", xmax(mesh))
        print(io, "    nelements: ", nelements(mesh))
    end
end

"""
    Mesh

Struct that holds the information for a simple homogeneous one-dimensional mesh.
"""
struct Mesh{NDIMS, RealT} <: AbstractMesh{NDIMS, RealT}
    xmin::RealT
    xmax::RealT
    N_elements::Int

    function Mesh{RealT}(xmin::RealT, xmax::RealT, N_elements::Int) where {RealT}
        @assert xmin < xmax
        @assert N_elements > 0
        new{1, RealT}(xmin, xmax, N_elements)
    end
end

"""
    Mesh(xmin, xmax, N_elements)

Create a simple homogeneous one-dimensional mesh from `xmin` to `xmax` with `N_elements` elements.
"""
function Mesh(xmin, xmax, N_elements)
    xmin, xmax = promote(xmin, xmax)
    Mesh{typeof(xmin)}(xmin, xmax, N_elements)
end

"""
    xmin(mesh::AbstractMesh)

Return the minimum coordinate of the mesh.
"""
xmin(mesh::Mesh) = mesh.xmin
"""
    xmax(mesh::AbstractMesh)

Return the maximum coordinate of the mesh.
"""
xmax(mesh::Mesh) = mesh.xmax
@inline nelements(mesh::Mesh) = mesh.N_elements

"""
    element_spacing(mesh::AbstractMesh)

Return the spacing of the elements in the mesh.
This is the length of each element, which is the same for all elements in a homogeneous [`Mesh`](@ref)
and can be different for each element in an [`InhomogeneousMesh`](@ref).
"""
element_spacing(mesh::Mesh) = (xmax(mesh) - xmin(mesh)) / nelements(mesh)

"""
    element_spacing(mesh::AbstractMesh, element)

Return the length of the element `element` in the mesh.
"""
element_spacing(mesh::AbstractMesh, element) = element_spacing(mesh)

"""
    left_element_boundary(mesh::AbstractMesh, element)

Return the left boundary coordinate of the element `element` in the mesh.
"""
function left_element_boundary(mesh::AbstractMesh, element)
    @assert 1<=element<=nelements(mesh) "Element index out of bounds"
    dx = element_spacing(mesh)
    return xmin(mesh) + (element - 1) * dx
end

"""
    InhomogeneousMesh{NDIMS, RealT}
    InhomogeneousMesh(coordinates)

Create an inhomogeneous one-dimensional mesh from the given `coordinates`.
"""
struct InhomogeneousMesh{NDIMS, RealT} <: AbstractMesh{NDIMS, RealT}
    coordinates::Vector{RealT}

    function InhomogeneousMesh(coordinates::AbstractVector{RealT}) where {RealT}
        @assert length(coordinates)>1 "InhomogeneousMesh requires at least two coordinates"
        new{1, RealT}(sort(collect(coordinates)))
    end
end

xmin(mesh::InhomogeneousMesh) = first(mesh.coordinates)
xmax(mesh::InhomogeneousMesh) = last(mesh.coordinates)
@inline nelements(mesh::InhomogeneousMesh) = length(mesh.coordinates) - 1
function element_spacing(mesh::InhomogeneousMesh, element)
    mesh.coordinates[element + 1] - mesh.coordinates[element]
end
function left_element_boundary(mesh::InhomogeneousMesh, element)
    @assert 1<=element<=nelements(mesh) "Element index out of bounds"
    return mesh.coordinates[element]
end

"""
    OversetGridMesh{NDIMS, RealT, MeshLeft, MeshRight}
    OversetGridMesh(mesh_left, mesh_right)

Create an overset grid (Chimera) mesh that combines two meshes, `mesh_left` and `mesh_right`, which have an overlap region.
"""
struct OversetGridMesh{NDIMS, RealT, MeshLeft <: AbstractMesh, MeshRight <: AbstractMesh} <:
       AbstractMesh{NDIMS, RealT}
    mesh_left::MeshLeft
    mesh_right::MeshRight

    function OversetGridMesh(mesh_left::MeshLeft,
                             mesh_right::MeshRight) where {RealT,
                                                           MeshLeft <:
                                                           AbstractMesh{1, RealT},
                                                           MeshRight <:
                                                           AbstractMesh{1, RealT}}
        @assert xmin(mesh_left) < xmin(mesh_right) &&
                xmax(mesh_left) > xmin(mesh_right) &&
                xmin(mesh_left) < xmax(mesh_right)
        new{1, RealT, MeshLeft, MeshRight}(mesh_left, mesh_right)
    end
end

xmin(mesh::OversetGridMesh) = xmin(mesh.mesh_left)
xmax(mesh::OversetGridMesh) = xmax(mesh.mesh_right)
@inline nelements(mesh::OversetGridMesh) = nelements(mesh.mesh_left) +
                                           nelements(mesh.mesh_right)
function element_spacing(mesh::OversetGridMesh, element)
    if element <= nelements(mesh.mesh_left)
        return element_spacing(mesh.mesh_left, element)
    else
        return element_spacing(mesh.mesh_right, element - nelements(mesh.mesh_left))
    end
end
function left_element_boundary(mesh::OversetGridMesh, element)
    if element <= nelements(mesh.mesh_left)
        return left_element_boundary(mesh.mesh_left, element)
    else
        return left_element_boundary(mesh.mesh_right, element - nelements(mesh.mesh_left))
    end
end

function left_overlap_element(mesh::OversetGridMesh)
    b = xmin(mesh.mesh_right)
    l_left = findfirst(element -> left_element_boundary(mesh.mesh_left, element) >= b,
                       1:nelements(mesh.mesh_left)) - 1
    return l_left
end
function right_overlap_element(mesh::OversetGridMesh)
    c = xmax(mesh.mesh_left)
    l_right = findfirst(element -> left_element_boundary(mesh.mesh_right, element) >= c,
                        1:nelements(mesh.mesh_right)) - 1
    return l_right
end

function Base.show(io::IO, ::MIME"text/plain",
                   mesh::OversetGridMesh{NDIMS, RealT}) where {NDIMS, RealT}
    if get(io, :compact, false)
        show(io, mesh)
    else
        println(io, string(nameof(typeof(mesh))), "{", RealT, "} ")
        println(io, "    mesh_left: ", mesh.mesh_left)
        print(io, "    mesh_right: ", mesh.mesh_right)
    end
end
