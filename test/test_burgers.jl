@testitem "burgers_basic.jl without source terms" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "burgers_basic.jl"),
                        source_terms=nothing,
                        l2=[0.757360668118935], linf=[0.8290216540696225],
                        cons_error=[3.228528555609955e-13],
                        change_mass=-3.228528555609955e-13,
                        change_entropy=-0.43782112131374884,
                        entropy_timederivative=-0.04501282859656586)
end

@testitem "burgers_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "burgers_basic.jl"),
                        l2=[0.004357967003300109], linf=[0.00933123559299709],
                        cons_error=[9.50350909079134e-14],
                        change_mass=-9.50350909079134e-14,
                        change_entropy=-0.00014813858633111465,
                        entropy_timederivative=-0.0002894716669563646)
end
