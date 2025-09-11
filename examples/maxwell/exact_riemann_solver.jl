using SimpleDiscontinuousGalerkin

equations = MaxwellEquations1D(0.5)

u_ll = SVector(-1.0, 2.0)
u_rr = SVector(2.0, 1.0)

prob = RiemannProblem(u_ll, u_rr)
solver = RiemannSolver(prob, equations)

t = 0.0:0.01:1.0
x = -1.0:0.01:1.0

sol = solve(solver, x, t)
