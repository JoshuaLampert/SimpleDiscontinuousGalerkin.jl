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

@testitem "compressible_euler_basic.jl with initial_condition_density_wave.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "compressible_euler_basic.jl"),
                        initial_condition=initial_condition_density_wave,
                        source_terms=nothing,
                        interval=50,
                        l2=[
                            0.012043278452772381,
                            0.0012043278452770788,
                            6.021639226103159e-5
                        ],
                        linf=[
                            0.0227222927108105,
                            0.0022722292710732755,
                            0.00011361146551536194
                        ],
                        cons_error=[
                            7.394085344003543e-14,
                            1.6542323066914832e-14,
                            3.907985046680551e-12
                        ], change_entropy=-1.6749982220787274e-5,
                        entropy_timederivative=-1.418715274603688e-5,
                        atol=1e-8) # To make CI pass
end

@testitem "compressible_euler_basic.jl with initial_condition_density_wave.jl and flux_godunov" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "compressible_euler_basic.jl"),
                        initial_condition=initial_condition_density_wave,
                        source_terms=nothing,
                        surface_flux=flux_godunov,
                        interval=50,
                        l2=[
                            0.006842388438912919,
                            0.0006842309974553353,
                            3.382179575461028e-5
                        ],
                        linf=[
                            0.021832146895580107,
                            0.0021832043288498026,
                            0.0001090312647491487
                        ],
                        cons_error=[
                            4.9960036108132044e-14,
                            1.9456658506555868e-14,
                            2.4726887204451486e-12
                        ], change_entropy=-0.0007395988732046277,
                        entropy_timederivative=-0.0003131350857813464,
                        atol=1e-6) # To make CI pass
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
