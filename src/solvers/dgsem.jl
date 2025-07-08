"""
    DGSEM(; RealT=Float64, polydeg::Integer,
            surface_flux=flux_central,
            surface_integral=SurfaceIntegralWeakForm(surface_flux),
            volume_integral=VolumeIntegralWeakForm())

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

@inline polydeg(solver::DGSEM) = length(grid(solver.basis)) - 1

Base.summary(io::IO, solver::DGSEM) = print(io, "DGSEM(polydeg=$(polydeg(solver)))")

"""
    FDSBP(D; RealT=Float64,
             surface_flux=flux_central,
             surface_integral=SurfaceIntegralWeakForm(surface_flux),
             volume_integral=VolumeIntegralWeakForm())

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

Base.summary(io::IO, solver::FDSBP) = print(io, "FDSBP(D=$(solver.basis)")

function create_jacobian_and_node_coordinates(mesh, solver::Union{DGSEM, FDSBP})
    nodes_basis = grid(solver.basis)
    dx_basis = last(nodes_basis) - first(nodes_basis) # length of the basis nodes
    # compute all mapped nodes
    x = zeros(real(solver), nnodes(solver), nelements(mesh))
    jacobian = zeros(real(solver), nelements(mesh))
    for element in eachelement(mesh)
        x_l = left_element_boundary(mesh, element)
        dx = element_spacing(mesh, element) # length of the element
        jacobian[element] = dx / dx_basis
        for j in eachindex(nodes_basis)
            x[j, element] = x_l + jacobian[element] * (nodes_basis[j] - first(nodes_basis))
        end
    end
    return jacobian, x
end

function create_tmp_scalar(mesh, solver)
    return VectorOfArray([zeros(real(solver), nnodes(solver, element))
                          for element in eachelement(mesh)])
end

function create_cache(mesh, equations, solver, initial_condition, boundary_conditions)
    jacobian, x = create_jacobian_and_node_coordinates(mesh, solver)
    tmp_scalar = create_tmp_scalar(mesh, solver)
    cache = (; jacobian, node_coordinates = x, tmp_scalar,
             create_cache(mesh, equations, solver, solver.volume_integral)...,
             create_cache(mesh, equations, solver, solver.surface_integral)...)
    return cache
end

function apply_jacobian!(du, mesh, equations, solver, cache)
    apply_jacobian!(du, mesh, equations, solver, cache, cache.jacobian)
end
function apply_jacobian!(du, mesh, equations, solver, cache, jacobian)
    for element in eachelement(mesh)
        for i in eachnode(solver, element)
            for v in eachvariable(equations)
                du[v, i, element] = du[v, i, element] / jacobian[element]
            end
        end
    end
end
