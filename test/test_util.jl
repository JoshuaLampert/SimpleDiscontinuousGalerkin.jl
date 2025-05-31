using Test: @test
using TrixiTest: @trixi_test_nowarn

# Use a macro to avoid world age issues when defining new initial conditions etc.
# inside an example.
"""
    @test_trixi_include(example; l2=nothing, linf = nothing,
                                 atol=1e-12, rtol=sqrt(eps()))

Test by calling `trixi_include(example; parameters...)`.
By default, only the absence of error output is checked.
If `l2` of `linf` are specified, in addition the resulting L2/Linf errors for each
variable are compared approximately against these reference values, using
`atol, rtol` as absolute/relative tolerance.
"""
macro test_trixi_include(example, args...)
    local l2 = get_kwarg(args, :l2, nothing)
    local linf = get_kwarg(args, :linf, nothing)
    local atol = get_kwarg(args, :atol, 1e-12)
    local rtol = get_kwarg(args, :rtol, sqrt(eps()))

    local kwargs = Pair{Symbol, Any}[]
    for arg in args
        if (arg.head == :(=) &&
            !(arg.args[1] in (:l2, :linf, :atol, :rtol)))
            push!(kwargs, Pair(arg.args...))
        end
    end

    quote
        println("═"^100)
        println($example)

        # evaluate examples in the scope of the module they're called from
        @trixi_test_nowarn trixi_include(@__MODULE__, $example; $kwargs...)

        if !isnothing($l2) || !isnothing($linf)
            # TODO: This should use the proper L2/Linf norms based on the basis of the solver, probably
            # also implemented as `AnalysisCallback`
            ini_cond = semi.initial_condition.(flat_grid(semi), sol.t[end], equations)
            for v in eachvariable(equations)
                diff = get_variable(sol.u[end], v, semi) .- getindex.(ini_cond, v)
                l2 = sqrt(sum(diff .^ 2) / (ndofs(semi)))
                linf = maximum(abs.(diff))
                if !isnothing($l2)
                    @test isapprox(l2, $l2[v], atol = $atol, rtol = $rtol)
                end
                if !isnothing($linf)
                    @test isapprox(linf, $linf[v], atol = $atol, rtol = $rtol)
                end
            end
        end
        println("═"^100)
    end
end

# Get the first value assigned to `keyword` in `args` and return `default_value`
# if there are no assignments to `keyword` in `args`.
function get_kwarg(args, keyword, default_value)
    val = default_value
    for arg in args
        if arg.head == :(=) && arg.args[1] == keyword
            val = arg.args[2]
            break
        end
    end
    return val
end
