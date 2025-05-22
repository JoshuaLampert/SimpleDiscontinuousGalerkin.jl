module SimpleDiscontinuousGalerkin

import LinearAlgebra: Diagonal, diag
using Reexport: @reexport
using SimpleUnPack: @unpack
@reexport using StaticArrays: SVector
@reexport using SummationByPartsOperators
import SummationByPartsOperators: AbstractDerivativeOperator, grid
using TimerOutputs: TimerOutputs, print_timer, reset_timer!
@reexport using TrixiBase: trixi_include
using TrixiBase: TrixiBase, @trixi_timeit, timer

include("equations/equations.jl")
include("mesh.jl")
include("boundary_conditions.jl")
include("solvers/solver.jl")
include("semidiscretization.jl")
include("callbacks_step/callbacks_step.jl")

export cons2cons
export LinearAdvectionEquation1D
export flux, flux_central, flux_godunov
export initial_condition_convergence_test
export Mesh
export boundary_condition_periodic, boundary_condition_do_nothing,
       BoundaryConditionDirichlet
export DGSEM, FDSBP,
       VolumeIntegralStrongForm, VolumeIntegralWeakForm,
       VolumeIntegralFluxDifferencing, VolumeIntegralFluxDifferencingStrongForm,
       SurfaceIntegralStrongForm, SurfaceIntegralWeakForm
export Semidiscretization, semidiscretize
export SummaryCallback
end
