@testitem "Maxwell_basic.jl" setup=[Setup] begin
    @test_trixi_include(joinpath(examples_dir(), "Maxwell_basic.jl"),
                        l2=[0.002295157400147163, 0.001538334261211901],
                        linf=[0.00491084212854303, 0.0026176927803420458],
                        cons_error=[1.6653345369377348e-16, 4.579669976578771e-16],
                        change_entropy=-3.473529797126673e-9,
                        entropy_timederivative=-0.18015669687179633)
end
