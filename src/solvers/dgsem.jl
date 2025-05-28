"""
    DGSEM(; RealT=Float64, polydeg::Integer,
            surface_flux=flux_central,
            surface_integral=SurfaceIntegralWeakForm(surface_flux),
            volume_integral=VolumeIntegralWeakForm(),
            mortar=MortarL2(basis))

Create a discontinuous Galerkin spectral element method (DGSEM) using a
`LegendreDerivativeOperator` with polynomials of degree `polydeg`.
"""
const DGSEM = DG{Basis} where {Basis <: LegendreDerivativeOperator}

# The constructor using only keyword arguments is convenient for elixirs since
# it allows to modify the polynomial degree and other parameters via
# `trixi_include`.
function DGSEM(; RealT = Float64,
               polydeg::Integer,
               surface_flux = flux_central,
               surface_integral = SurfaceIntegralWeakForm(surface_flux),
               volume_integral = VolumeIntegralWeakForm())
    N = polydeg + 1
    xmin = RealT(-1.0)
    xmax = RealT(1.0)
    basis = legendre_derivative_operator(xmin, xmax, N)
    return DG{typeof(basis), typeof(surface_integral),
              typeof(volume_integral)}(basis, surface_integral, volume_integral)
end

@inline polydeg(dg::DGSEM) = length(grid(dg.basis)) - 1

Base.summary(io::IO, dg::DGSEM) = print(io, "DGSEM(polydeg=$(polydeg(dg)))")

"""
    FDSBP(D; RealT=Float64,
             surface_flux=flux_central,
             surface_integral=SurfaceIntegralWeakForm(surface_flux),
             volume_integral=VolumeIntegralWeakForm(),
             mortar=MortarL2(basis))

Create a discontinuous Galerkin method using a summation-by-parts operator
`D` from SummationByPartsOperators.jl. This is similar to the `DGSEM`, but uses
a general derivative operator instead of a Legendre derivative operator.
"""
const FDSBP = DG{Basis} where {Basis <: AbstractDerivativeOperator}

function FDSBP(D; RealT = Float64,
               surface_flux = flux_central,
               surface_integral = SurfaceIntegralWeakForm(surface_flux),
               volume_integral = VolumeIntegralWeakForm())
    basis = D
    return DG{typeof(basis), typeof(surface_integral),
              typeof(volume_integral)}(basis, surface_integral, volume_integral)
end

Base.summary(io::IO, dg::FDSBP) = print(io, "FDSBP(D=$D)")

function create_cache(mesh, equations, dg::Union{DGSEM, FDSBP}, initial_condition,
                      boundary_conditions)
    dx = element_spacing(mesh) # length of each element
    nodes_basis = grid(dg.basis)
    dx_basis = last(nodes_basis) - first(nodes_basis) # length of the basis nodes
    jacobian = dx / dx_basis
    # compute all mapped nodes
    x = zeros(real(dg), nnodes(dg), nelements(mesh))
    for element in eachelement(mesh)
        x_l = xmin(mesh) + (element - 1) * dx
        for j in eachindex(nodes_basis)
            x[j, element] = x_l + jacobian * (nodes_basis[j] - first(nodes_basis))
        end
    end
    cache = (; jacobian, node_coordinates = x,
             create_cache(mesh, equations, dg, dg.volume_integral)...,
             create_cache(mesh, equations, dg, dg.surface_integral)...)
    return cache
end

function apply_jacobian!(du, mesh, equations, dg::Union{DGSEM, FDSBP}, cache)
    (; jacobian) = cache
    @. du = du / jacobian
end
