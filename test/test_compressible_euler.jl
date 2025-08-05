@testitem "compressible_euler_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "compressible_euler_basic.jl"),
                        l2=[
                            2.9956346259676573e-5,
                            1.8489334163606683e-5,
                            7.110148970433093e-5
                        ],
                        linf=[
                            6.125900116904504e-5,
                            5.023528949710254e-5,
                            0.00016930866446784876
                        ],
                        cons_error=[
                            2.708944180085382e-14,
                            2.353672812205332e-14,
                            4.796163466380676e-14
                        ], change_entropy=3.5522082164618496e-8,
                        entropy_timederivative=3.304866663320083e-8)
end
