module SimpleDiscontinuousGalerkin

import LinearAlgebra: Diagonal
using Reexport: @reexport
using SimpleUnPack: @unpack
@reexport using StaticArrays: SVector
@reexport using SummationByPartsOperators
import SummationByPartsOperators: grid
@reexport using TrixiBase: trixi_include
using TrixiBase: TrixiBase, @trixi_timeit, timer

include("equations/equations.jl")
include("mesh.jl")
include("boundary_conditions.jl")
include("solvers/solver.jl")
include("semidiscretization.jl")

export LinearAdvectionEquation1D, flux_central, flux_godunov,
       Mesh, grid, boundary_condition_periodic,
       DGSEM, Semidiscretization, semidiscretize
end
