"""
    SimpleDiscontinuousGalerkin

**SimpleDiscontinuousGalerkin.jl** is a Julia package that implements some basic discontinuous Galerkin (DG) methods for the solution of
hyperbolic partial differential equations (PDEs). The package is designed to be simple and easy to use and understand. It is intended for
educational purposes and to provide a starting point for more complex DG methods. For a more comprehensive and advanced implementation of
DG methods, we recommend using the package [Trixi.jl](https://github.com/trixi-framework/Trixi.jl). This package can be understood as a
minimalistic version of Trixi.jl, which is designed to be easy to understand and modify. Many design concepts are inspired by Trixi.jl,
but the implementation is much simpler and more straightforward.

SimpleDiscontinuousGalerkin.jl builds on the foundations of the package [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).

See also: [SimpleDiscontinuousGalerkin.jl](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl)
"""
module SimpleDiscontinuousGalerkin

import LinearAlgebra: Diagonal, diag
using RecursiveArrayTools: VectorOfArray
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

export cons2cons, eachvariable, nvariables
export LinearAdvectionEquation1D
export flux, flux_central, flux_godunov
export initial_condition_convergence_test
export Mesh, nelements, eachelement, element_spacing
export boundary_condition_periodic, boundary_condition_do_nothing,
       BoundaryConditionDirichlet
export eachnode, nnodes, nelements, ndofs, flat_grid, get_variable
export DGSEM, FDSBP, PerElementFDSBP,
       VolumeIntegralStrongForm, VolumeIntegralWeakForm,
       VolumeIntegralFluxDifferencing, VolumeIntegralFluxDifferencingStrongForm,
       SurfaceIntegralStrongForm, SurfaceIntegralWeakForm
export Semidiscretization, semidiscretize
export SummaryCallback
end
