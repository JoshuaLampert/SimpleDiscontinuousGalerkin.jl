@testitem "Aqua.jl" begin
    import Aqua
    using ExplicitImports: test_explicit_imports

    # Aqua.jl
    Aqua.test_all(SimpleDiscontinuousGalerkin)

    # ExplicitImports.jl
    # Non-public names from SummationByPartsOperators.jl, which we need to dispatch
    # on the operator types
    non_public_explicit_imports = (:AbstractDerivativeOperator,)
    non_public_qualified_accesses = (
                                     # SummationByPartsOperators.jl
                                     :AbstractNonperiodicDerivativeOperator, :xmin, :xmax,
                                     # PolynomialBases.jl
                                     :barycentric_weights,
                                     # Base, used in the `AnalysisCallback`
                                     :gc_live_bytes, :gc_time_ns,
                                     # TimerOutputs.jl
                                     :complement!,
                                     # TrixiBase.jl
                                     :find_assignment)
    # Names can only be declared `public` since Julia v1.11. Thus, we only check
    # whether the imported and accessed names are public on Julia v1.11 and newer.
    check_public = VERSION >= v"1.11"
    imports_are_public = check_public && (ignore = non_public_explicit_imports,)
    accesses_are_public = check_public && (ignore = non_public_qualified_accesses,)
    test_explicit_imports(SimpleDiscontinuousGalerkin;
                          # We reexport SummationByPartsOperators.jl
                          no_implicit_imports = (;
                                                 skip = (Core, Base,
                                                         SimpleDiscontinuousGalerkin.SummationByPartsOperators)),
                          all_explicit_imports_are_public = imports_are_public,
                          all_qualified_accesses_are_public = accesses_are_public)
end
