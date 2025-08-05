@testitem "maxwell_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "maxwell_basic.jl"),
                        l2=[0.0014845551649291684, 0.0018613613260282855],
                        linf=[0.003311775580174814, 0.0029746676945732464],
                        cons_error=[2.498001805406602e-16, 5.342948306008566e-16],
                        change_entropy=-0.000112697123739125,
                        entropy_timederivative=-0.00011456197705196625)
end
