@testitem "maxwell_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_MAXWELL, "maxwell_basic.jl"),
                        l2=[0.0010497390241669513, 0.0013161812158726952],
                        linf=[0.00331177558017548, 0.002974667694574537],
                        cons_error=[2.498001805406602e-16, 5.342948306008566e-16],
                        change_entropy=-0.000112697123739125,
                        entropy_timederivative=-0.00011456197705196625)
end
