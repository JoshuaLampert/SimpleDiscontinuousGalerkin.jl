using Documenter
using SimpleDiscontinuousGalerkin
using TrixiBase

# Dynamically set all files in subdirectories of the source directory to include all files in these subdirectories
# This way they don't need to be listed explicitly
EQUATIONS_FILES = joinpath.(Ref("equations"),
                            readdir(joinpath(dirname(@__DIR__), "src",
                                             "equations")))
CALLBACKS_STEP_FILES = joinpath.(Ref("callbacks_step"),
                                 readdir(joinpath(dirname(@__DIR__), "src",
                                                  "callbacks_step")))
SOLVER_FILES = joinpath.(Ref("solvers"),
                          readdir(joinpath(dirname(@__DIR__), "src",
                                           "solvers")))

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(SimpleDiscontinuousGalerkin, :DocTestSetup,
                    :(using SimpleDiscontinuousGalerkin);
                    recursive = true)
DocMeta.setdocmeta!(TrixiBase, :DocTestSetup, :(using TrixiBase); recursive = true)

makedocs(;
         modules = [SimpleDiscontinuousGalerkin, TrixiBase],
         authors = "Joshua Lampert <joshua.lampert@uni-hamburg.de>",
         repo = Remotes.GitHub("JoshuaLampert", "SimpleDiscontinuousGalerkin.jl"),
         sitename = "SimpleDiscontinuousGalerkin.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://JoshuaLampert.github.io/SimpleDiscontinuousGalerkin.jl/stable",
                                  edit_link = "main"),
         pages = ["Home" => "index.md",
             "Development" => "development.md",
             "Reference" => [
                 "TrixiBase" => "ref-trixibase.md",
                 "SimpleDiscontinuousGalerkin" => "ref.md"
             ],
             "License" => "license.md"])

deploydocs(;
           repo = "github.com/JoshuaLampert/SimpleDiscontinuousGalerkin.jl",
           devbranch = "main",
           push_preview = true)
