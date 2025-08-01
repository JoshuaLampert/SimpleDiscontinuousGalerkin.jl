# Development

If you have any suggestions or ideas for improvements or new features, we are pleased to accept and discuss
[issues](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/issues) or if you are willing to contribute,
feel free to [open a pull request](https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl/pulls), even if it
is only fixing a typo or improving the docs.

## Changing SimpleDiscontinuousGalerkin.jl and running it locally

If you plan to edit SimpleDiscontinuousGalerkin.jl, you first need to clone a local copy of the repository, which can
be done by using `git`. It is recommended that you create a project, e.g. call it `run`, inside the repository,
where you can add packages that you use during executing and testing SimpleDiscontinuousGalerkin.jl, but are not needed
by SimpleDiscontinuousGalerkin.jl. This way you can keep the Project.toml of the main repository clean. To do so, you
can execute the following lines in a terminal:

```sh
git clone https://github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl.git
cd SimpleDiscontinuousGalerkin
mkdir run
cd run
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path=".."))' # Install local SimpleDiscontinuousGalerkin.jl clone
julia --project=. -e 'using Pkg; Pkg.add(["OrdinaryDiffEqLowStorageRK", "Plots"])' # Install additional packages
```

If you use other packages for executing SimpleDiscontinuousGalerkin.jl, you can add them to the project in the `run`
directory in an analogous way as above. To use the Julia project within `run`, be sure to start the Julia REPL
by

```sh
julia --project=.
```

if already inside the the `run` directory or `julia --project=run` if in the main directory of the repo.

## Preview of the documentation

If you want to build the documentation locally, you can run

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

once from the SimpleDiscontinuousGalerkin.jl main directory to tell [Documenter.jl](https://documenter.juliadocs.org/stable/man/guide/)
to build the documentation of your local clone. To build the documentation, run

```sh
julia --project=docs --color=yes docs/make.jl
```

The resulting `.html` files can then be found in `docs/build/` and you can look at them by opening them in a browser.
For pull requests from the main repository (i.e. not from a fork), the documentation is automatically built and can
be previewed under `https://JoshuaLampert.github.io/SimpleDiscontinuousGalerkin.jl/previews/PRXXX/` where `XXX` is the number
of the pull request.

## Testing

In the test suite of SimpleDiscontinuousGalerkin.jl, the results of the examples are tested against reference solutions.
This means that if you add an example and verified it works as expected, you should put the reference solution
into the corresponding test. To obtain the reference solution and automatically copy it into the clipboard, you can run

```julia
errs = errors(analysis_callback)
ints = integrals(analysis_callback)
using InteractiveUtils
clipboard("l2=$(errs.l2_error[:, end])\n, linf=$(errs.linf_error[:, end])\n, cons_error=$(errs.conservation_error[:, end])\n, change_mass=$(ints.mass[end] - ints.mass[1]),\nchange_entropy=$(ints.entropy[end] - ints.entropy[1]), \nentropy_timederivative=$(ints.entropy_timederivative[end])")
```

Note that some of the quantities are only computed if they are included in the `analysis_callback`, which means that you
might need to adjust the above command depending on which quantities you want to test for.
