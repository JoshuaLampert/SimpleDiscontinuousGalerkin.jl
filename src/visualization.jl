struct PlotData{Conversion}
    semisol::Pair{<:Semidiscretization, <:ODESolution}
    plot_initial::Bool
    conversion::Conversion
    step::Integer
end

@recipe function f(semi::Semidiscretization, data, data_exact, names, plot_initial)
    x = flat_grid(semi)
    for i in 1:length(names)
        if plot_initial == true
            @series begin
                subplot --> i
                label --> names[i]
                xguide --> "x"
                yguide --> names[i]
                title --> names[i]
                x, data_exact[i, :]
            end
        end
        @series begin
            subplot --> i
            linestyle --> :solid
            label --> "initial $(names[i])"
            x, data[i, :]
        end
    end
end

@recipe function f(semi::SemidiscretizationOversetGrid, data, data_exact, names,
                   plot_initial)
    x_left, x_right = flat_grid(semi)
    for i in 1:length(names)
        if plot_initial == true
            data_exact_left, data_exact_right = data_exact
            @series begin
                subplot --> i
                linestyle --> :solid
                label --> "initial $(names[i])"
                x_left, data_exact_left[i, :]
            end
            @series begin
                subplot --> i
                linestyle --> :solid
                label --> "initial $(names[i])"
                x_right, data_exact_right[i, :]
            end
        end
        data_left, data_right = data
        @series begin
            subplot --> i
            label --> names[i]
            xguide --> "x"
            yguide --> names[i]
            title --> names[i]
            x_left, data_left[i, :]
        end
        @series begin
            subplot --> i
            label --> names[i]
            xguide --> "x"
            yguide --> names[i]
            title --> names[i]
            x_right, data_right[i, :]
        end
    end
end

function compute_data(mesh::OversetGridMesh, equations, solver, u, u_exact, plot_initial,
                      conversion, nvars)
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    u_left, u_right = u
    u_exact_left, u_exact_right = u_exact

    data_left, data_left_exact = compute_data(mesh_left, equations, solver_left, u_left,
                                              u_exact_left,
                                              plot_initial, conversion, nvars)
    data_right, data_right_exact = compute_data(mesh_right, equations, solver_right,
                                                u_right, u_exact_right,
                                                plot_initial, conversion, nvars)
    return (data_left, data_right), (data_left_exact, data_right_exact)
end

function compute_data(mesh, equations, solver, u, u_exact, plot_initial, conversion, nvars)
    data = zeros(real(mesh), nvars, ndofs(mesh, solver))
    if plot_initial == true
        data_exact = zeros(real(mesh), nvars, ndofs(mesh, solver))
    else
        data_exact = nothing
    end
    j = 1
    for element in eachelement(mesh)
        for node in eachnode(solver, element)
            data[:, j] .= conversion(get_node_vars(u, equations, node, element), equations)
            if plot_initial == true
                data_exact[:, j] .= conversion(get_node_vars(u_exact, equations, node,
                                                             element), equations)
            end
            j += 1
        end
    end
    return data, data_exact
end

@recipe function f(plotdata::PlotData)
    @unpack semisol, plot_initial, conversion, step = plotdata
    semi, sol = semisol
    mesh, equations, solver, _ = mesh_equations_solver_cache(semi)
    names = varnames(conversion, equations)

    nvars = length(conversion(zeros(nvariables(semi)), equations))

    if step == -1
        step = length(sol.t)
    end

    initial_condition = semi.initial_condition
    t = sol.t[step]
    u = sol.u[step]
    if plot_initial
        u_exact = compute_coefficients(initial_condition, t, semi)
    else
        u_exact = nothing
    end

    data, data_exact = compute_data(mesh, equations, solver, u, u_exact, plot_initial,
                                    conversion, nvars)

    plot_title --> "$(get_name(equations)) at t = $(round(t, digits=5))"
    layout --> nvars

    semi, data, data_exact, names, plot_initial
end

@recipe function f(semisol::Pair{<:Semidiscretization, <:ODESolution}; plot_initial = false,
                   conversion = cons2prim, step = -1)
    PlotData(semisol, plot_initial, conversion, step)
end

@recipe function f(semi::Semidiscretization, sol::ODESolution; plot_initial = false,
                   conversion = cons2prim, step = -1)
    PlotData(semi => sol, plot_initial, conversion, step)
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

@recipe function f(riemann_solution::RiemannSolverSolution; step = -1,
                   conversion = cons2cons)
    equations = riemann_solution.solver.equations
    names = varnames(conversion, equations)

    nvars = length(conversion(zeros(nvariables(equations)), equations))
    nsubplots = nvars

    if step == -1
        step = length(riemann_solution.solution)
    end

    t = riemann_solution.t[step]
    x = riemann_solution.x
    data = conversion.(riemann_solution[step], equations)
    plot_title --> "$(get_name(equations)) at t = $(round(t, digits=5))"
    layout --> nsubplots
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
