using SimpleDiscontinuousGalerkin

equations = BurgersEquation1D()

u_ll = -1.0
u_rr = 2.0

prob = RiemannProblem(u_ll, u_rr)
solver = RiemannSolver(prob, equations)

t = 0.0:0.01:0.5
x = -1.0:0.01:1.0

sol = solve(solver, x, t)
