module SimpleDiscontinuousGalerkin

import LinearAlgebra: Diagonal
using Reexport: @reexport
using SimpleUnPack: @unpack
@reexport using StaticArrays: SVector
@reexport using SummationByPartsOperators
import SummationByPartsOperators: AbstractDerivativeOperator, grid
@reexport using TrixiBase: trixi_include
using TrixiBase: TrixiBase, @trixi_timeit, timer

include("equations/equations.jl")
include("mesh.jl")
include("boundary_conditions.jl")
include("solvers/solver.jl")
include("semidiscretization.jl")

export LinearAdvectionEquation1D
export flux_central, flux_godunov
export initial_condition_convergence_test
export Mesh, grid
export boundary_condition_periodic, boundary_condition_do_nothing,
       BoundaryConditionDirichlet
export DGSEM, FDSBP, VolumeIntegralStrongForm, VolumeIntegralWeakForm,
       SurfaceIntegralStrongForm, SurfaceIntegralWeakForm
export Semidiscretization, semidiscretize
end
