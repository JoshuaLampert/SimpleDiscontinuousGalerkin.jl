SemidiscretizationOversetGrid = Semidiscretization{<:OversetGridMesh}

digest_solver(::OversetGridMesh, solver::Union{DGSEM, FDSBP}) = (solver, solver)

function ndofs(equations, mesh::OversetGridMesh, solver::Tuple)
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    return ndofs(equations, mesh_left, solver_left) +
           ndofs(equations, mesh_right, solver_right)
end

function allocate_coefficients(mesh::OversetGridMesh, equations, solver::Tuple)
    solver_left, solver_right = solver
    u_left = allocate_coefficients(mesh.mesh_left, equations, solver_left)
    u_right = allocate_coefficients(mesh.mesh_right, equations, solver_right)
    return VectorOfArray([u_left, u_right])
end

function compute_coefficients!(u, func, t, mesh::OversetGridMesh, equations, solver, cache)
    u_left, u_right = u
    solver_left, solver_right = solver
    cache_left, cache_right = cache
    compute_coefficients!(u_left, func, t, mesh.mesh_left, equations, solver_left,
                          cache_left)
    compute_coefficients!(u_right, func, t, mesh.mesh_right, equations, solver_right,
                          cache_right)
end

function grid(semi::SemidiscretizationOversetGrid)
    cache_left, cache_right = semi.cache
    return (cache_left.node_coordinates, cache_right.node_coordinates)
end

function flat_grid(semi::SemidiscretizationOversetGrid)
    return collect(Iterators.flatten(parent(grid(semi)[1]))),
           collect(Iterators.flatten(parent(grid(semi)[2])))
end

function get_variable(u, v, semi::SemidiscretizationOversetGrid)
    u_left, u_right = u
    solver_left, solver_right = semi.solver
    return get_variable(u_left, v, solver_left), get_variable(u_right, v, solver_right)
end

function Iterators.flatten(semi::SemidiscretizationOversetGrid, u::VectorOfArray)
    solver_left, solver_right = semi.solver
    u_left, u_right = u
    u_left_flattened = Iterators.flatten(solver_left, u_left)
    u_right_flattened = Iterators.flatten(solver_right, u_right)
    return [collect(u_left_flattened); collect(u_right_flattened)]
end

# function jacobian_fd(semi::SemidiscretizationOversetGrid;
#                      t0 = zero(real(semi)),
#                      u0_ode = compute_coefficients(semi.initial_condition, t0, semi))
#     # copy the initial state since it will be modified in the following
#     u_ode = copy(u0_ode)
#     du0_ode = similar(u_ode)
#     dup_ode = similar(u_ode)
#     dum_ode = similar(u_ode)

#     # compute residual of linearization state
#     rhs!(du0_ode, u_ode, semi, t0)

#     # initialize Jacobian matrix
#     J = zeros(eltype(u_ode), ndofs(semi), ndofs(semi))

#     idx = 0
#     # use second order finite difference to estimate Jacobian matrix
#     meshes = [semi.mesh.mesh_left, semi.mesh.mesh_right]
#     for nmesh in 1:2
#         for element in eachelement(meshes[nmesh])
#             for node in eachnode(semi.solver[nmesh], element)
#                 for v in eachvariable(semi.equations)
#                     idx += 1
#                     # determine size of fluctuation
#                     epsilon = sqrt(eps(typeof(u0_ode[v, node, element, nmesh])))

#                     # plus fluctuation
#                     u_ode[v, node, element, nmesh] = u0_ode[v, node, element, nmesh] +
#                                                      epsilon
#                     rhs!(dup_ode, u_ode, semi, t0)

#                     # minus fluctuation
#                     u_ode[v, node, element, nmesh] = u0_ode[v, node, element, nmesh] -
#                                                      epsilon
#                     rhs!(dum_ode, u_ode, semi, t0)

#                     # restore linearization state
#                     u_ode[v, node, element, nmesh] = u0_ode[v, node, element, nmesh]

#                     # central second order finite difference
#                     Ji_ode = (dup_ode .- dum_ode) ./ (2 * epsilon)
#                     Ji = flatten(semi, Ji_ode)
#                     @. J[:, idx] = Ji
#                 end
#             end
#         end
#     end
#     return J
# end

function create_cache(mesh::OversetGridMesh, equations, solver)
    solver_left, solver_right = solver
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    cache_left = create_cache(mesh_left, equations, solver_left)
    cache_right = create_cache(mesh_right, equations, solver_right)
    return (; cache_left, cache_right)
end

function rhs!(du, u, t, mesh::OversetGridMesh, equations, initial_condition,
              boundary_conditions, source_terms, solver, cache)
    du_left, du_right = du
    u_left, u_right = u
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    cache_left, cache_right = cache

    @trixi_timeit timer() "reset ∂u/∂t" begin
        reset_du!(du_left)
        reset_du!(du_right)
    end

    @trixi_timeit timer() "volume integral" begin
        calc_volume_integral!(du_left, u_left, mesh_left, equations,
                              solver_left.volume_integral, solver_left, cache_left)
        calc_volume_integral!(du_right, u_right, mesh_right, equations,
                              solver_right.volume_integral, solver_right, cache_right)
    end

    @trixi_timeit timer() "interface flux" begin
        calc_interface_flux!(cache_left.surface_flux_values, u_left, mesh_left,
                             equations, solver_left.surface_integral, solver_left,
                             cache_left)
        calc_interface_flux!(cache_right.surface_flux_values, u_right, mesh_right,
                             equations, solver_right.surface_integral, solver_right,
                             cache_right)
    end

    @trixi_timeit timer() "boundary flux" begin
        surface_flux_values = (cache_left.surface_flux_values,
                               cache_right.surface_flux_values)
        surface_integral = (solver_left.surface_integral, solver_right.surface_integral)
        calc_boundary_flux!(surface_flux_values, u, t,
                            boundary_conditions, mesh, equations,
                            surface_integral, solver)
    end

    @trixi_timeit timer() "surface integral" begin
        calc_surface_integral!(du_left, u_left, mesh_left, equations,
                               solver_left.surface_integral, solver_left, cache_left)
        calc_surface_integral!(du_right, u_right, mesh_right, equations,
                               solver_right.surface_integral, solver_right, cache_right)
    end

    @trixi_timeit timer() "Jacobian" begin
        apply_jacobian!(du_left, mesh_left, equations, solver_left, cache_left)
        apply_jacobian!(du_right, mesh_right, equations, solver_right, cache_right)
    end

    @trixi_timeit timer() "source terms" begin
        calc_sources!(du_left, u_left, t, source_terms, mesh_left, equations,
                      solver_left, cache_left)
        calc_sources!(du_right, u_right, t, source_terms, mesh_right, equations,
                      solver_right, cache_right)
    end
end

function calc_boundary_flux!(surface_flux_values, u, t, boundary_conditions,
                             mesh::OversetGridMesh, equations,
                             integral::Tuple, solver)
    surface_flux_values_left, surface_flux_values_right = surface_flux_values
    u_left, u_right = u
    (; x_neg, x_pos) = boundary_conditions
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    integral_left, integral_right = integral
    solver_left, solver_right = solver

    # Left boundary condition of left mesh
    u_ll = x_neg(u, xmin(mesh), t, mesh, equations, solver, true)
    u_rr = get_node_vars(u_left, equations, 1, 1)
    f = integral_left.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values_left, f, equations, 1, 1)

    # Right boundary condition of left mesh
    u_ll = get_node_vars(u_left, equations, nnodes(solver_left, nelements(mesh_left)),
                         nelements(mesh_left))
    # TODO: Make this more efficient by computing e_M_right and e_M_left only once
    # during initialization and storing it in the cache
    l_right = right_overlap_element(mesh)
    c = xmax(mesh_left)
    D = get_basis(solver_right, l_right)
    linear_map(x, a, b, c, d) = c + (x - a) / (b - a) * (d - c)
    xlr_L = left_element_boundary(mesh_right, l_right)
    xlr_R = left_element_boundary(mesh_right, l_right + 1)
    c_mapped = linear_map(c, xlr_L, xlr_R, first(grid(D)), last(grid(D)))
    u_rr = zeros(real(solver_right), nvariables(equations))
    e_M_right = interpolation_operator(c_mapped, D)
    for v in eachvariable(equations)
        values = u_right[v, :, l_right]
        u_rr[v] = e_M_right' * values
    end
    f = integral_left.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values_left, f, equations, 2, nelements(mesh_left))

    # Left boundary condition of right mesh
    l_left = left_overlap_element(mesh)
    b = xmin(mesh_right)
    D = get_basis(solver_left, l_left)
    xll_L = left_element_boundary(mesh_left, l_left)
    xll_R = left_element_boundary(mesh_left, l_left + 1)
    b_mapped = linear_map(b, xll_L, xll_R, first(grid(D)), last(grid(D)))
    u_ll = zeros(real(solver_left), nvariables(equations))
    e_M_left = interpolation_operator(b_mapped, D)
    for v in eachvariable(equations)
        values = u_left[v, :, l_left]
        u_ll[v] = e_M_left' * values
    end
    u_rr = get_node_vars(u_right, equations, 1, 1)
    f = integral_right.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values_right, f, equations, 1, 1)

    # Right boundary condition of right mesh
    u_ll = get_node_vars(u_right, equations, nnodes(solver_right, nelements(mesh_right)),
                         nelements(mesh_right))
    u_rr = x_pos(u, xmax(mesh), t, mesh, equations, solver, false)
    f = integral_right.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values_right, f, equations, 2, nelements(mesh_right))
end

# This method is for integrating a vector quantity for all variables over the entire domain,
# such as the whole solution vector `u` (`VectorOfArray{T, 4, Vector{Array{T, 3}}}` for DG methods
# with same basis across elements and `VectorOfArray{T, 4, Vector{VectorOfArray{T, 3, Vector{Matrix{T}}}}}}`
# for `PerElementFDSBP`).
function PolynomialBases.integrate(func,
                                   u::VectorOfArray{T, 4, T1},
                                   semi::SemidiscretizationOversetGrid) where {T, T1}
    u_left, u_right = u
    integrals = zeros(real(semi), nvariables(semi))
    mesh_left, mesh_right = semi.mesh.mesh_left, semi.mesh.mesh_right
    for v in eachvariable(semi)
        u_left_v = VectorOfArray([u_left[v, :, element]
                                  for element in eachelement(mesh_left)])
        u_right_v = VectorOfArray([u_right[v, :, element]
                                   for element in eachelement(mesh_right)])
        u_v = VectorOfArray([u_left_v, u_right_v])
        integrals[v] = integrate(func, u_v, semi)
    end
    return integrals
end

# This method is for integrating a scalar quantity over the entire domain.
# Need to dispatch on type of `u` to avoid method ambiguities.
function PolynomialBases.integrate(func,
                                   u::VectorOfArray{T, 3, T1},
                                   semi::SemidiscretizationOversetGrid) where {T,
                                                                               T1 <:
                                                                               Vector{VectorOfArray{T,
                                                                                                    2,
                                                                                                    Vector{Vector{T}}}}
                                                                               }
    u_left, u_right = u
    jacobian_left, jacobian_right = semi.cache.cache_left.jacobian,
                                    semi.cache.cache_right.jacobian
    mesh_left, mesh_right = semi.mesh.mesh_left, semi.mesh.mesh_right
    solver_left, solver_right = semi.solver
    l_left = left_overlap_element(semi.mesh)
    # TODO: This integrates the left domain only until the right boundary of the left overlap element,
    # i.e., we count the integral from b to the right boundary of the left overlap element twice.
    # TODO: Which part should be integrated in the overlap region? This uses the right mesh,
    # which is not correct for negative velocities.
    integral_left = sum(integrate_on_element(func, u_left.u[element],
                                             get_basis(solver_left, element), element,
                                             jacobian_left)
                        for element in 1:l_left)
    integral_right = sum(integrate_on_element(func, u_right.u[element],
                                              get_basis(solver_right, element),
                                              element, jacobian_right)
                         for element in 1:nelements(mesh_right))
    return integral_left + integral_right
end

function integrate_quantity(func, u, semi::SemidiscretizationOversetGrid)
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    u_left, u_right = u
    cache_left, cache_right = cache
    quantity_left, quantity_right = get_tmp_cache_scalar(cache_left),
                                    get_tmp_cache_scalar(cache_right)
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    compute_quantity!(quantity_left, func, u_left, mesh_left, equations, solver_left)
    compute_quantity!(quantity_right, func, u_right, mesh_right, equations, solver_right)
    quantity = VectorOfArray([quantity_left, quantity_right])
    integrate(quantity, semi)
end

function analyze(::typeof(entropy_timederivative), du, u, t,
                 semi::SemidiscretizationOversetGrid)
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    u_left, u_right = u
    du_left, du_right = du
    cache_left, cache_right = cache
    quantity_left, quantity_right = get_tmp_cache_scalar(cache_left),
                                    get_tmp_cache_scalar(cache_right)
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    compute_quantity_timederivative!(quantity_left, cons2entropy, du_left, u_left,
                                     mesh_left, equations, solver_left)
    compute_quantity_timederivative!(quantity_right, cons2entropy, du_right, u_right,
                                     mesh_right, equations, solver_right)
    quantity = VectorOfArray([quantity_left, quantity_right])
    return integrate(quantity, semi)
end

function calc_error_norms(u, t, initial_condition, mesh::OversetGridMesh,
                          equations, solver, cache)
    u_left, u_right = u
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    cache_left, cache_right = cache
    l2_error_left, linf_error_left = calc_error_norms(u_left, t, initial_condition,
                                                      mesh_left, equations,
                                                      solver_left, cache_left)
    l2_error_right, linf_error_right = calc_error_norms(u_right, t, initial_condition,
                                                        mesh_right, equations,
                                                        solver_right, cache_right)
    return l2_error_left + l2_error_right, max(linf_error_left, linf_error_right)
end

function max_dt(u, t, mesh::OversetGridMesh, equations, solver, cache)
    u_left, u_right = u
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    cache_left, cache_right = cache
    max_scaled_speed_left = max_dt(u_left, t, mesh_left, equations, solver_left,
                                   cache_left)
    max_scaled_speed_right = max_dt(u_right, t, mesh_right, equations, solver_right,
                                    cache_right)
    return max(max_scaled_speed_left, max_scaled_speed_right)
end
