struct PlotData
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    plot_initial::Bool
    step::Integer
end

@recipe function f(plotdata::PlotData)
    @unpack semisol, plot_initial, step = plotdata
    semi, sol = semisol
    equations = semi.equations
    names = varnames(cons2cons, equations)

    nvars = nvariables(semi)
    nsubplots = nvars

    if step == -1
        step = length(sol.t)
    end

    initial_condition = semi.initial_condition
    t = sol.t[step]

    if plot_initial == true
        u_exact = compute_coefficients(initial_condition, t, semi)
    end

    u = sol.u[step]

    plot_title --> "$(get_name(equations)) at t = $(round(t, digits=5))"
    layout --> nsubplots

    for i in 1:nsubplots
        if plot_initial == true
            if semi isa SemidiscretizationOversetGrid
                u_exact_left, u_exact_right = get_variable(u_exact, i, semi)
                x_left, x_right = flat_grid(semi)
                @series begin
                    subplot --> i
                    linestyle --> :solid
                    label --> "initial $(names[i])"
                    x_left, u_exact_left
                end
                @series begin
                    subplot --> i
                    linestyle --> :solid
                    label --> "initial $(names[i])"
                    x_right, u_exact_right
                end
            else
                @series begin
                    subplot --> i
                    linestyle --> :solid
                    label --> "initial $(names[i])"
                    flat_grid(semi), get_variable(u_exact, i, semi)
                end
            end
        end

        if semi isa SemidiscretizationOversetGrid
            u_left, u_right = get_variable(u, i, semi)
            x_left, x_right = flat_grid(semi)
            @series begin
                subplot --> i
                label --> names[i]
                xguide --> "x"
                yguide --> names[i]
                title --> names[i]
                x_left, u_left
            end
            @series begin
                subplot --> i
                label --> names[i]
                xguide --> "x"
                yguide --> names[i]
                title --> names[i]
                x_right, u_right
            end
        else
            @series begin
                subplot --> i
                label --> names[i]
                xguide --> "x"
                yguide --> names[i]
                title --> names[i]
                flat_grid(semi), get_variable(u, i, semi)
            end
        end
    end
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}; plot_initial = false,
                   step = -1)
    PlotData(semisol, plot_initial, step)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution; plot_initial = false,
                   step = -1)
    PlotData(semi => sol, plot_initial, step)
end

function pretty_form_utf(name)
    if name == :l2_error
        return "L² error"
    elseif name == :linf_error
        return "L∞ error"
    elseif name == :conservation_error
        return "∫|u_u₀|"
    else
        return string(name)
    end
end

@recipe function f(cb::DiscreteCallback{Condition, Affect!}; what = (:integrals,),
                   label_extension = "", start_from = 1,
                   exclude = (:entropy_timederivative,)) where {Condition,
                                                                Affect! <: AnalysisCallback}
    t = tstops(cb)
    @assert length(t)>start_from "The keyword argument `start_from` needs to be smaller than the number of timesteps: $(length(t))"
    subplot = 1
    layout --> length(what)
    if :integrals in what
        ints = integrals(cb)

        for (i, (name, integral)) in enumerate(pairs(ints))
            name in exclude && continue
            quantity = cb.affect!.analysis_integrals[i]
            @series begin
                subplot --> subplot
                label --> pretty_form_utf(quantity) * " " * label_extension
                title --> "change of invariants"
                xguide --> "t"
                yguide --> "change of invariants"
                t[start_from:end], (integral .- integral[1])[start_from:end]
            end
        end
        subplot += 1
    end
    if :errors in what
        errs = errors(cb)
        for (name, err) in pairs(errs)
            name in exclude && continue
            @series begin
                subplot --> subplot
                label --> pretty_form_utf(name) * " " * label_extension
                title --> "errors"
                xguide --> "t"
                yguide --> "sum of errors"
                t[start_from:end], dropdims(sum(err, dims = 1), dims = 1)[start_from:end]
            end
        end
    end
end

@recipe function f(riemann_solution::RiemannSolverSolution; step = -1)
    equations = riemann_solution.solver.equations
    names = varnames(cons2cons, equations)

    nvars = nvariables(equations)
    nsubplots = nvars

    if step == -1
        step = length(riemann_solution.solution)
    end

    t = riemann_solution.t[step]
    x = riemann_solution.x
    data = riemann_solution[step]
    plot_title --> "$(get_name(equations)) at t = $(round(t, digits=5))"
    for i in 1:nsubplots
        @series begin
            subplot --> i
            label --> names[i]
            xguide --> "x"
            yguide --> names[i]
            title --> names[i]
            x, getindex.(data, i)
        end
    end
end