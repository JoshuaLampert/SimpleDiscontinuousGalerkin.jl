using SimpleDiscontinuousGalerkin

equations = LinearAdvectionEquation1D(-2.0)

u_ll = SVector(1.0)
u_rr = SVector(0.0)

prob = RiemannProblem(u_ll, u_rr)
solver = RiemannSolver(equations)

t = 0.0:0.1:1.0
x = -1.0:0.1:1.0

sol = solve(prob, solver, x, t)
