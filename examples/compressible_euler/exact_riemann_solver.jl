using SimpleDiscontinuousGalerkin

equations = CompressibleEulerEquations1D(1.4)

# Sod's shock tube
u_ll = prim2cons(SVector(1.0, 0.0, 1.0), equations)
u_rr = prim2cons(SVector(0.125, 0.0, 0.1), equations)

prob = RiemannProblem(u_ll, u_rr)
solver = RiemannSolver(equations)

t = 0.0:0.01:0.5
x = -1.0:0.01:1.0

sol = solve(prob, solver, x, t)
