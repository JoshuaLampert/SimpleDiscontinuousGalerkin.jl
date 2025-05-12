# SimpleDiscontinuousGalerkin.jl

[![Build Status](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/graph/badge.svg?token=ZnS5D3tWSK)](https://codecov.io/gh/JoshuaLampert/SimpleDiscontinuousGalerkin.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

**SimpleDiscontinuousGalerkin.jl** is a [Julia](https://julialang.org/) package that
implements some basic discontinuous Galerkin (DG) methods for the solution of hyperbolic
partial differential equations (PDEs). The package is designed to be simple and easy to use and understand. It is
intended for educational purposes and to provide a starting point for more complex DG methods. For a more
comprehensive and advanced implementation of DG methods, we recommend using the package
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl). This package can be understood as a
minimalistic version of Trixi.jl, which is designed to be easy to understand and modify. Many design concepts are
inspired by Trixi.jl, but the implementation is much simpler and more straightforward. SimpleDiscontinuousGalerkin.jl
builds on the foundations of the package [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl).

## Installation

If you have not yet installed Julia, then you first need to [download Julia](https://julialang.org/downloads/). Please
[follow the instructions for your operating system](https://julialang.org/downloads/platform/). SimpleDiscontinuousGalerkin.jl
works with Julia v1.10 and newer. You can install SimpleDiscontinuousGalerkin.jl by
executing the following commands from the Julia REPL

```julia
julia> using Pkg

julia> Pkg.add("https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl")
```

## Usage

In the Julia REPL, first load the package SimpleDiscontinuousGalerkin.jl

```julia
julia> using SimpleDiscontinuousGalerkin
```

SimpleDiscontinuousGalerkin.jl is built on top of the package SummationByPartsOperators.jl and exports all the functions
and types of the package.

## Authors

The package is developed and maintained by Joshua Lampert (University of Hamburg).

## License and contributing

SimpleDiscontinuousGalerkin.jl is published under the MIT license (see [License](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/blob/main/LICENSE)).
We are pleased to accept contributions from everyone, preferably in the form of a PR.
