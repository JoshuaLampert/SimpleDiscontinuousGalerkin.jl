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
                   exclude = (:entropy_timederivative,)) where {Condition, Affect! <: AnalysisCallback}
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
                label := pretty_form_utf(quantity) * " " * label_extension
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
