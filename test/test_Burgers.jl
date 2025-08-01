@testitem "Burgers_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "Burgers_basic.jl"),
                        l2=[0.757360668118935], linf=[0.8290216540696225],
                        cons_error=[3.228528555609955e-13],
                        change_mass=-3.228528555609955e-13,
                        change_entropy=-0.43782112131374884,
                        entropy_timederivative=-0.04501282859656586)
end
