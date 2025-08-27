function vandermonde_matrix(functions, nodes)
    N = length(nodes)
    K = length(functions)
    T = typeof(functions[1](nodes[1]))
    V = zeros(T, N, K)
    for i in 1:N
        for j in 1:K
            V[i, j] = functions[j](nodes[i])
        end
    end
    return V
end
function basis_functions(D::SummationByPartsOperators.AbstractNonperiodicDerivativeOperator)
    p = accuracy_order(D)
    return [x -> x^i for i in 0:p]
end
function lagrange_polynomial(xs, i)
    return x -> prod((x - xs[j]) / (xs[i] - xs[j]) for j in 1:length(xs) if j != i)
end
function basis_functions(D::LegendreDerivativeOperator)
    return [lagrange_polynomial(grid(D), i) for i in 1:length(grid(D))]
end
function interpolation_operator(x,
                                D::SummationByPartsOperators.AbstractNonperiodicDerivativeOperator)
    nodes = grid(D)
    functions = basis_functions(D)
    V = vandermonde_matrix(functions, nodes)
    values = [functions[i](x) for i in eachindex(functions)]
    e_M = V' \ values
    return e_M
end
# Collocation (e.g., Gauss-Lobatto)
function interpolate(x, values,
                     D::SummationByPartsOperators.AbstractNonperiodicDerivativeOperator)
    e_M = interpolation_operator(x, D)
    return e_M' * values
end
# No collocation (e.g., Finite Difference)
function interpolate(x, values, D::SummationByPartsOperators.DerivativeOperator)
    nodes = grid(D)
    spl = Spline1D(nodes, values)
    return spl(x)
end

@inline function ln_mean(x::RealT, y::RealT) where {RealT <: Real}
    epsilon_f2 = convert(RealT, 1.0e-4)
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        return (x + y) / @evalpoly(f2,
                         2,
                         convert(RealT, 2 / 3),
                         convert(RealT, 2 / 5),
                         convert(RealT, 2 / 7))
    else
        return (y - x) / log(y / x)
    end
end

@inline function inv_ln_mean(x::RealT, y::RealT) where {RealT <: Real}
    epsilon_f2 = convert(RealT, 1.0e-4)
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        return @evalpoly(f2,
                         2,
                         convert(RealT, 2 / 3),
                         convert(RealT, 2 / 5),
                         convert(RealT, 2 / 7)) / (x + y)
    else
        return log(y / x) / (y - x)
    end
end

"""
    examples_dir()

Return the directory where the example files provided by SimpleDiscontinuousGalerkin.jl are located. If SimpleDiscontinuousGalerkin.jl is
installed as a regular package (with `]add SimpleDiscontinuousGalerkin`), these files are read-only and should *not* be
modified. To find out which files are available, use, e.g., `readdir`.

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Examples
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(SimpleDiscontinuousGalerkin, "examples")::String

"""
    default_example()

Return the path to an example that can be used to quickly see SimpleDiscontinuousGalerkin.jl in action.
See also [`examples_dir`](@ref).

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
"""
function default_example()
    joinpath(examples_dir(), "linear_advection", "linear_advection.jl")
end

function convergence_test(example, iterations_or_Ns; kwargs...)
    convergence_test(Main, example, iterations_or_Ns; kwargs...)
end

"""
    convergence_test([mod::Module=Main,] example, iterations; io::IO = stdout, kwargs...)
    convergence_test([mod::Module=Main,] example, Ns::AbstractVector; io::IO = stdout, kwargs...)

Run multiple simulations using the setup given in `example` and compute
the experimental order of convergence (EOC) in the ``L^2`` and ``L^\\infty`` norm.
If `iterations` is passed as integer, in each iteration, the resolution of the respective mesh
will be doubled. If `Ns` is passed as vector, the simulations will be run for each value of `Ns`.
Additional keyword arguments `kwargs...` and the optional module `mod` are passed directly
to [`trixi_include`](@ref).

Adjusted from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
"""
function convergence_test(mod::Module, example, iterations; io::IO = stdout,
                          kwargs...)
    @assert(iterations>1,
            "Number of iterations must be bigger than 1 for a convergence analysis")

    initial_N = extract_initial_N(example, kwargs)
    Ns = initial_N * 2 .^ (0:(iterations - 1))
    convergence_test(mod, example, Ns; io = io, kwargs...)
end

function convergence_test(mod::Module, example, Ns::AbstractVector; io::IO = stdout,
                          kwargs...)
    # Types of errors to be calculated
    errors = Dict(:l2 => Float64[], :linf => Float64[])

    Base.require_one_based_indexing(Ns)
    sort!(Ns)
    iterations = length(Ns)
    # run simulations and extract errors
    for iter in 1:iterations
        println("Running convtest iteration ", iter, "/", iterations)

        trixi_include(mod, example; kwargs..., N_elements = Ns[iter])

        l2_error, linf_error = mod.analysis_callback(mod.sol)

        # collect errors as one vector to reshape later
        append!(errors[:l2], l2_error)
        append!(errors[:linf], linf_error)

        println("\n\n")
        println("#"^100)
    end

    # Use raw error values to compute EOC
    analyze_convergence(io, errors, iterations, mod.semi, Ns)
end

# Analyze convergence for any semidiscretization
# Note: this intermediate method is to allow dispatching on the semidiscretization
function analyze_convergence(io, errors, iterations, semi, Ns)
    _, equations, _, _ = mesh_equations_solver_cache(semi)
    variablenames = varnames(cons2cons, equations)
    analyze_convergence(io, errors, iterations, variablenames, Ns)
end

# This method is called with the collected error values to actually compute and print the EOC
function analyze_convergence(io, errors, iterations,
                             variablenames::Union{Tuple, AbstractArray}, Ns)
    nvariables = length(variablenames)

    # Reshape errors to get a matrix where the i-th row represents the i-th iteration
    # and the j-th column represents the j-th variable
    errorsmatrix = Dict(kind => transpose(reshape(error, (nvariables, iterations)))
                        for (kind, error) in errors)

    # Calculate EOCs where the columns represent the variables
    eocs = Dict(kind => log.(error[2:end, :] ./ error[1:(end - 1), :]) ./
                        log.(Ns[1:(end - 1)] ./ Ns[2:end])
                for (kind, error) in errorsmatrix)

    eoc_mean_values = Dict{Symbol, Any}()
    eoc_mean_values[:variables] = variablenames

    for (kind, error) in errorsmatrix
        println(io, kind)

        for v in variablenames
            @printf(io, "%-25s", v)
        end
        println(io, "")

        for k in 1:nvariables
            @printf(io, "%-5s", "N")
            @printf(io, "%-10s", "error")
            @printf(io, "%-10s", "EOC")
        end
        println(io, "")

        # Print errors for the first iteration
        for k in 1:nvariables
            @printf(io, "%-5d", Ns[1])
            @printf(io, "%-10.2e", error[1, k])
            @printf(io, "%-10s", "-")
        end
        println(io, "")

        # For the following iterations print errors and EOCs
        for j in 2:iterations
            for k in 1:nvariables
                @printf(io, "%-5d", Ns[j])
                @printf(io, "%-10.2e", error[j, k])
                @printf(io, "%-10.2f", eocs[kind][j - 1, k])
            end
            println(io, "")
        end
        println(io, "")

        # Print mean EOCs
        mean_values = zeros(nvariables)
        for v in 1:nvariables
            mean_values[v] = sum(eocs[kind][:, v]) ./ length(eocs[kind][:, v])
            @printf(io, "%-15s", "mean")
            @printf(io, "%-10.2f", mean_values[v])
        end
        eoc_mean_values[kind] = mean_values
        println(io, "")
        println(io, "-"^100)
    end

    return eoc_mean_values, errorsmatrix
end

function extract_initial_N(example, kwargs)
    code = read(example, String)
    expr = Meta.parse("begin \n$code \nend")

    if haskey(kwargs, :N_elements)
        return kwargs[:N_elements]
    else
        # get N from the example
        N = TrixiBase.find_assignment(expr, :N_elements)
        return N
    end
end
