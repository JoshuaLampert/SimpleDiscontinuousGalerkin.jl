@testitem "compressible_euler_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_EULER, "compressible_euler_basic.jl"),
                        l2=[
                            2.1182335579730997e-5,
                            1.3073933566041213e-5,
                            5.027634552313215e-5
                        ],
                        linf=[
                            6.125900117659455e-5,
                            5.0235289478894884e-5,
                            0.00016930866444653248
                        ],
                        cons_error=[
                            2.708944180085382e-14,
                            2.353672812205332e-14,
                            4.796163466380676e-14
                        ], change_entropy=3.5522082164618496e-8,
                        entropy_timederivative=3.304866663320083e-8)
end

@testitem "compressible_euler_basic.jl with flux_hll" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_EULER, "compressible_euler_basic.jl"),
                        surface_flux=flux_hll,
                        l2=[
                            1.4335766272237236e-5,
                            1.0668807342292978e-5,
                            4.1528318426135346e-5
                        ],
                        linf=[
                            4.764287292591263e-5,
                            4.505781308372647e-5,
                            0.00015040034055857632
                        ],
                        cons_error=[
                            2.708944180085382e-14,
                            2.708944180085382e-14,
                            5.3290705182007514e-14
                        ], change_entropy=3.407981097325319e-8,
                        entropy_timederivative=3.162079266830209e-8)
end

@testitem "compressible_euler_basic.jl with initial_condition_density_wave.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_EULER, "compressible_euler_basic.jl"),
                        initial_condition=initial_condition_density_wave,
                        source_terms=nothing,
                        interval=50,
                        l2=[
                            0.008515883851203871,
                            0.0008515883851205352,
                            4.257941924627521e-5
                        ],
                        linf=[
                            0.022722292756255147,
                            0.0022722292756178625,
                            0.00011361146577115733
                        ],
                        cons_error=[
                            7.394085344003543e-14,
                            1.6542323066914832e-14,
                            3.907985046680551e-12
                        ], change_entropy=-1.6749982220787274e-5,
                        entropy_timederivative=-1.418715274603688e-5,
                        atol=1e-8) # To make CI pass
end

@testitem "compressible_euler_basic.jl with initial_condition_density_wave.jl and flux_hll" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_EULER, "compressible_euler_basic.jl"),
                        surface_flux=flux_hll,
                        initial_condition=initial_condition_density_wave,
                        source_terms=nothing,
                        interval=50,
                        l2=[
                            0.008511707905745933,
                            0.0008511707905748125,
                            4.2558539526967456e-5
                        ],
                        linf=[
                            0.022718034169821655,
                            0.0022718034169656037,
                            0.00011359017280199168
                        ],
                        cons_error=[
                            7.527312106958561e-14,
                            1.6542323066914832e-14,
                            3.765876499528531e-12
                        ], change_entropy=-1.696940963391569e-5,
                        entropy_timederivative=-1.4441696482736521e-5,
                        atol=1e-8) # To make CI pass
end

@testitem "compressible_euler_basic.jl with initial_condition_density_wave.jl and flux_godunov" setup=[Setup] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR_EULER, "compressible_euler_basic.jl"),
                        initial_condition=initial_condition_density_wave,
                        source_terms=nothing,
                        surface_flux=flux_godunov,
                        interval=50,
                        l2=[
                            0.004838299661374554,
                            0.00048383400177452295,
                            2.424169595757841e-5
                        ],
                        linf=[
                            0.021832144910775897,
                            0.002183233771653803,
                            0.00010912357123515903
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
    @test_trixi_include(joinpath(EXAMPLES_DIR_EULER, "compressible_euler_basic.jl"),
                        initial_condition=initial_condition_weak_blast_wave,
                        source_terms=nothing,
                        volume_flux=flux_kennedy_gruber,
                        l2=[0.10872397221262474, 0.26939372211118634, 0.4054025830648324],
                        linf=[0.1861764323552595, 0.4910781570555981, 0.692006422257144],
                        cons_error=[
                            7.460698725481052e-14,
                            7.7021722333370235e-16,
                            1.8829382497642655e-13
                        ], change_entropy=-0.003335960418462871,
                        entropy_timederivative=-0.0029203971291180525)
end
