@testitem "Aqua.jl" begin
    import Aqua
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports,
                           check_all_explicit_imports_via_owners,
                           check_all_qualified_accesses_via_owners,
                           check_no_self_qualified_accesses

    Aqua.test_all(SimpleDiscontinuousGalerkin)
    @test isnothing(check_no_implicit_imports(SimpleDiscontinuousGalerkin,
                                              skip = (Core, Base,
                                                      SimpleDiscontinuousGalerkin.SummationByPartsOperators)))
    @test isnothing(check_no_stale_explicit_imports(SimpleDiscontinuousGalerkin))
    @test isnothing(check_all_qualified_accesses_via_owners(SimpleDiscontinuousGalerkin;
                                                            ignore = (:ustrip,)))
    @test isnothing(check_no_self_qualified_accesses(SimpleDiscontinuousGalerkin))
end
