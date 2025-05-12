"""
    Mesh

Struct that holds the information for a simple homogeneous one-dimensional mesh.
"""
struct Mesh{NDIMS, RealT}
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

function Base.show(io::IO, mesh::Mesh{RealT}) where {RealT}
    print(io, "Mesh{", RealT, "} with ", mesh.N_elements, " elements")
    print(io, " from ", mesh.xmin, " to ", mesh.xmax)
end

function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh{RealT}) where {RealT}
    if get(io, :compact, false)
        show(io, mesh)
    else
        println(io, "Mesh{", RealT, "} ")
        println(io, "    xmin: ", mesh.xmin)
        println(io, "    xmax: ", mesh.xmax)
        print(io, "    N_elements: ", mesh.N_elements)
    end
end

@inline Base.ndims(::Mesh{NDIMS}) where {NDIMS} = NDIMS
@inline nelements(mesh::Mesh) = mesh.N_elements
@inline eachelement(mesh::Mesh) = Base.OneTo(nelements(mesh))
@inline Base.real(::Mesh{NDIMS, RealT}) where {NDIMS, RealT} = RealT

xmin(mesh::Mesh) = mesh.xmin
xmax(mesh::Mesh) = mesh.xmax
