using Test: @test
using TrixiTest: @trixi_test_nowarn

# Use a macro to avoid world age issues when defining new initial conditions etc.
# inside an example.
"""
    @test_trixi_include(example; l2=nothing, linf = nothing, cons_error=nothing,
                        change_mass=nothing, change_entropy=nothing,
                        entropy_timederivative=nothing,
                        atol=1e-12, rtol=sqrt(eps()))

Test by calling `trixi_include(example; parameters...)`.
By default, only the absence of error output is checked.
If `l2`, `linf`, or `cons_error` are specified, in addition the resulting L2/Linf errors for each
variable are compared approximately against these reference values, using
`atol, rtol` as absolute/relative tolerance. If `change_mass` or `change_entropy` are specified,
the change in mass or entropy is compared against the reference values. If `entropy_timederivative`
is specified, the time derivative of the entropy is compared against the reference value
"""
macro test_trixi_include(example, args...)
    local l2 = get_kwarg(args, :l2, nothing)
    local linf = get_kwarg(args, :linf, nothing)
    local cons_error = get_kwarg(args, :cons_error, nothing)
    local change_mass = get_kwarg(args, :change_mass, nothing)
    local change_entropy = get_kwarg(args, :change_entropy, nothing)
    local entropy_timederivative = get_kwarg(args, :entropy_timederivative, nothing)
    local atol = get_kwarg(args, :atol, 1e-12)
    local rtol = get_kwarg(args, :rtol, sqrt(eps()))

    local kwargs = Pair{Symbol, Any}[]
    for arg in args
        if (arg.head == :(=) &&
            !(arg.args[1] in (:l2, :linf, :cons_error, :change_mass, :change_entropy,
                              :entropy_timederivative, :atol, :rtol)))
            push!(kwargs, Pair(arg.args...))
        end
    end

    quote
        println("═"^100)
        println($example)

        # evaluate examples in the scope of the module they're called from
        @trixi_test_nowarn trixi_include(@__MODULE__, $example; $kwargs...)

        # if present, compare l2, linf and conservation errors against reference values
        if !isnothing($l2) || !isnothing($linf) || !isnothing($cons_error)
            errs = errors(analysis_callback)

            if !isnothing($l2)
                l2_measured = errs.l2_error[:, end]
                @test length($l2) == length(l2_measured)
                for (l2_expected, l2_actual) in zip($l2, l2_measured)
                    @test isapprox(l2_expected, l2_actual, atol = $atol, rtol = $rtol)
                end
            end

            if !isnothing($linf)
                linf_measured = errs.linf_error[:, end]
                @test length($linf) == length(linf_measured)
                for (linf_expected, linf_actual) in zip($linf, linf_measured)
                    @test isapprox(linf_expected, linf_actual, atol = $atol, rtol = $rtol)
                end
            end

            if !isnothing($cons_error)
                cons_error_measured = errs.conservation_error[:, end]
                @test length($cons_error) == length(cons_error_measured)
                for (conservation_error_expected, conservation_error_actual) in zip($cons_error,
                                                                                    cons_error_measured)
                    @test isapprox(conservation_error_expected, conservation_error_actual,
                                   atol = $atol, rtol = $rtol)
                end
            end
        end

        if !isnothing($change_mass) || !isnothing($change_entropy) ||
           !isnothing($entropy_timederivative)
            ints = integrals(analysis_callback)

            if !isnothing($change_mass)
                mass_change_measured = ints.mass[end] - ints.mass[1]
                @test isapprox($change_mass, mass_change_measured,
                               atol = $atol, rtol = $rtol)
            end

            if !isnothing($change_entropy)
                entropy_change_measured = ints.entropy[end] - ints.entropy[1]
                @test isapprox($change_entropy, entropy_change_measured,
                               atol = $atol, rtol = $rtol)
            end

            if !isnothing($entropy_timederivative)
                entropy_timederivative_measured = ints.entropy_timederivative[end]
                @test isapprox($entropy_timederivative, entropy_timederivative_measured,
                               atol = $atol, rtol = $rtol)
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
