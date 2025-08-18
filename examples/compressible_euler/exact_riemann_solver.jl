using SimpleDiscontinuousGalerkin

equations = CompressibleEulerEquations1D(1.4)

# Sod's shock tube (t_end = 0.5)
u_ll = prim2cons(SVector(1.0, 0.0, 1.0), equations)
u_rr = prim2cons(SVector(0.125, 0.0, 0.1), equations)
# 123 problem (t_end = 0.3)
# u_ll = prim2cons(SVector(1.0, -2.0, 0.4), equations)
# u_rr = prim2cons(SVector(1.0, 2.0, 0.4), equations)
# Colliding flows (t_end = 0.4)
# u_ll = prim2cons(SVector(1.0, 3.0, 1.0), equations)
# u_rr = prim2cons(SVector(1.0, -3.0, 1.0), equations)
# Woodward-Colella (left part) (t_end = 0.024)
# u_ll = prim2cons(SVector(1.0, 0.0, 1000.0), equations)
# u_rr = prim2cons(SVector(1.0, 0.0, 0.01), equations)
# Woodward-Colella (left part) (t_end = 0.07)
# u_ll = prim2cons(SVector(1.0, 0.0, 0.01), equations)
# u_rr = prim2cons(SVector(1.0, 0.0, 100.0), equations)
# Shock collision (t_end = 0.07)
# u_ll = prim2cons(SVector(5.99924, 19.5975, 460.894), equations)
# u_rr = prim2cons(SVector(5.99242, -6.19633, 46.095), equations)

prob = RiemannProblem(u_ll, u_rr)
solver = RiemannSolver(prob, equations)

t = 0.0:0.01:0.5
x = -1.0:0.01:1.0

sol = solve(solver, x, t)
