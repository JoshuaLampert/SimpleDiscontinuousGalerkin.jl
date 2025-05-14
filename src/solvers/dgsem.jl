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
    dx = (xmax(mesh) - xmin(mesh)) / nelements(mesh) # length of each element
    # compute all mapped GLL nodes
    x = zeros(real(dg), nnodes(dg), nelements(mesh))
    for element in eachelement(mesh)
        x_l = xmin(mesh) + (element - 1) * dx + dx / 2
        for (j, xi_j) in enumerate(grid(dg))
            x[j, element] = x_l + dx / 2 * xi_j
        end
    end
    cache = (; dx, node_coordinates = x,
             create_cache(mesh, equations, dg, dg.volume_integral)...,
             create_cache(mesh, equations, dg, dg.surface_integral)...)
    return cache
end

function apply_jacobian!(du, mesh, equations, dg::Union{DGSEM, FDSBP}, cache)
    (; dx) = cache
    @. du = (2 / dx) * du
end
