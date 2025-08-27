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

using CommonSolve: CommonSolve, init, solve!
using Dierckx: Spline1D
import LinearAlgebra: Diagonal, diag, dot
using PolynomialBases: PolynomialBases
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe, @series
using RecursiveArrayTools: VectorOfArray
using Reexport: @reexport
using Roots: find_zero, AlefeldPotraShi
import SciMLBase: ODESolution, solve, u_modified!, get_tmp_cache, set_proposed_dt!
using SimpleUnPack: @unpack
@reexport using StaticArrays: SVector
@reexport using SummationByPartsOperators
import SummationByPartsOperators: AbstractDerivativeOperator, grid
using TimerOutputs: TimerOutputs, @notimeit, print_timer, reset_timer!
@reexport using TrixiBase: trixi_include
using TrixiBase: TrixiBase, @trixi_timeit, timer

include("utils.jl")
include("equations/equations.jl")
include("mesh.jl")
include("boundary_conditions.jl")
include("solvers/solver.jl")
include("semidiscretization/semidiscretization.jl")
include("callbacks_step/callbacks_step.jl")
include("visualization.jl")

export examples_dir, default_example, convergence_test
export cons2cons, cons2entropy, cons2prim, prim2cons, eachvariable, nvariables, varnames
export mass, entropy, electric_field, magnetic_field, density, velocity, momentum, pressure,
       density_pressure, entropy_thermodynamic, entropy_math, energy_total, energy_kinetic,
       energy_internal, energy_internal_specific
export LinearAdvectionEquation1D, BurgersEquation1D, MaxwellEquations1D,
       CompressibleEulerEquations1D
export FluxLaxFriedrichs, FluxHLL, flux, flux_central, flux_godunov, flux_lax_friedrichs,
       flux_ec, flux_hll, flux_ranocha, flux_kennedy_gruber
export RiemannProblem, RiemannSolver
export initial_condition_convergence_test, source_terms_convergence_test,
       initial_condition_density_wave, initial_condition_weak_blast_wave
export Mesh, InhomogeneousMesh, OversetGridMesh, nelements, eachelement, element_spacing
export boundary_condition_periodic, boundary_condition_do_nothing,
       BoundaryConditionDirichlet
export eachnode, nnodes, nelements, ndofs, flat_grid, get_variable
export DGSEM, FDSBP, PerElementFDSBP,
       VolumeIntegralStrongForm, VolumeIntegralWeakForm,
       VolumeIntegralFluxDifferencing, VolumeIntegralFluxDifferencingStrongForm,
       SurfaceIntegralStrongForm, SurfaceIntegralWeakForm
export Semidiscretization, semidiscretize, jacobian_fd
export SummaryCallback, AnalysisCallback, StepsizeCallback
export tstops, errors, integrals
export init, solve!, solve
end
