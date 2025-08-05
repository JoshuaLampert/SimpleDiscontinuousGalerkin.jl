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

@testitem "compressible_euler_basic.jl with initial_condition_weak_blast_wave.jl and flux_kennedy_gruber" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "compressible_euler_basic.jl"),
                        initial_condition=initial_condition_weak_blast_wave,
                        source_terms=nothing,
                        volume_flux=flux_kennedy_gruber,
                        l2=[0.15375891605816916, 0.3809802554278062, 0.5733258311913755],
                        linf=[0.18617643235524484, 0.4910781570555968, 0.6920064222571098],
                        cons_error=[
                            7.460698725481052e-14,
                            7.7021722333370235e-16,
                            1.8829382497642655e-13
                        ], change_entropy=-0.003335960418462871,
                        entropy_timederivative=-0.0029203971291180525)
end
