# Overview

This package provides a simple implementation of Discontinuous Galerkin (DG) methods for solving hyperbolic
conservation laws. It is designed to be easy to understand and extend, making it suitable for educational purposes
and quick prototyping of new ideas. It is not intended for high-performance computing or production use.
Many ideas and concepts are inspired by the package [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) and
most users will likely want to use that package instead. For the basic concepts of Discontinuous Galerkin methods,
which are mostly the same as in Trixi.jl, please refer to the [Trixi.jl documentation](https://trixi-framework.github.io/TrixiDocumentation/stable/).
The syntax of the API is very similar to Trixi.jl. Please check out
[the examples/ folder](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/tree/main/examples)
for some usage examples.

## Differences to Trixi.jl

- The focus of SimpleDiscontinuousGalerkin.jl is on understandability and extendability of the code and less on
  performance and features. It is designed closely follow the mathematical description of the DG method, which can
  make it suitable for teaching and learning purposes. For example, SimpleDiscontinuousGalerkin.jl implements more
  equivalent formulations like weak and strong formulation of the flux-differencing method to demonstrate the
  equivalence of these formulations.
- SimpleDiscontinuousGalerkin.jl has no support for many features Trixi.jl has (many more equations, multiple space
  dimensions, more sophisticated numerical fluxes, adaptive mesh refinement, shock capturing, positive-presreving
  methods, limiting strategies, callbacks, parallel computing, etc.).
- SimpleDiscontinuousGalerkin.jl has some (rather niche) features Trixi.jl does not support (yet), which stem from
  quick tests and prototyping to test some new ideas. This includes different numerical fluxes at the boundary than
  in the interior and different bases for each element via the [`PerElementBasis`](@ref) type.
- In contrast to Trixi.jl, which implements an own Gauss-Lobatto-Legendre basis for the Discontinuous Galerkin spectral
  element method ([`DGSEM`](@ref)), this package highlights more the underlying summation-by-parts (SBP) structure of
  the DG method. It reuses the `legendre_derivative_operator` from
  [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl) for the Discontinuous
  Galerkin spectral element method ([`DGSEM`](@ref)). This can be seen as a special case of [`FDSBP`](@ref) with
  the underlying SBP operator being the Legendre derivative operator.
