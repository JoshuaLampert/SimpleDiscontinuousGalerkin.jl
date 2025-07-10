function vandermonde_matrix(functions, nodes)
    N = length(nodes)
    K = length(functions)
    T = typeof(functions[1](nodes[1]))
    V = zeros(T, N, K)
    for i in 1:N
        for j in 1:K
            V[i, j] = functions[j](nodes[i])
        end
    end
    return V
end
function basis_functions(D::SummationByPartsOperators.AbstractNonperiodicDerivativeOperator)
    p = accuracy_order(D)
    return [x -> x^i for i in 0:p]
end
function lagrange_polynomial(xs, i)
    return x -> prod((x - xs[j]) / (xs[i] - xs[j]) for j in 1:length(xs) if j != i)
end
function basis_functions(D::LegendreDerivativeOperator)
    return [lagrange_polynomial(grid(D), i) for i in 1:length(grid(D))]
end
function interpolation_operator(x,
                                D::SummationByPartsOperators.AbstractNonperiodicDerivativeOperator)
    nodes = grid(D)
    functions = basis_functions(D)
    V = vandermonde_matrix(functions, nodes)
    values = [functions[i](x) for i in eachindex(functions)]
    e_M = V' \ values
    return e_M
end
