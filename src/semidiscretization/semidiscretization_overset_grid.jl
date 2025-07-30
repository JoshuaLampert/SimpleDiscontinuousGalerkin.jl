SemidiscretizationOversetGrid = Semidiscretization{<:OversetGridMesh}

digest_solver(::OversetGridMesh, solver::Union{DGSEM, FDSBP}) = (solver, solver)

# This allows us to treat a `Tuple` of `Solver`s as a `Solver`.
function Base.getproperty(solver::Tuple{S1, S2}, name::Symbol) where {S1 <: DG, S2 <: DG}
    return getproperty(solver[1], name)
end

function ndofs(mesh::OversetGridMesh, solver::Tuple)
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    return ndofs(mesh_left, solver_left) + ndofs(mesh_right, solver_right)
end

function compute_integral_operator(mesh, solver::Tuple, integral; kwargs...)
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    return (compute_integral_operator(mesh_left, solver_left, integral; kwargs...),
            compute_integral_operator(mesh_right, solver_right, integral; kwargs...))
end

function create_jacobian_and_node_coordinates(mesh::OversetGridMesh, solver::Tuple)
    solver_left, solver_right = solver
    jacobian_left, x_left = create_jacobian_and_node_coordinates(mesh.mesh_left,
                                                                 solver_left)
    jacobian_right, x_right = create_jacobian_and_node_coordinates(mesh.mesh_right,
                                                                   solver_right)
    return (jacobian_left, jacobian_right), (x_left, x_right)
end

function create_tmp_scalar(mesh::OversetGridMesh, solver)
    solver_left, solver_right = solver
    tmp_scalar_left = create_tmp_scalar(mesh.mesh_left, solver_left)
    tmp_scalar_right = create_tmp_scalar(mesh.mesh_right, solver_right)
    return VectorOfArray([tmp_scalar_left, tmp_scalar_right])
end

function allocate_coefficients(mesh::OversetGridMesh, equations, solver::Tuple)
    solver_left, solver_right = solver
    u_left = allocate_coefficients(mesh.mesh_left, equations, solver_left)
    u_right = allocate_coefficients(mesh.mesh_right, equations, solver_right)
    return VectorOfArray([u_left, u_right])
end

function compute_coefficients!(u, func, t, mesh::OversetGridMesh, equations, solver,
                               cache, node_coordinates)
    u_left, u_right = u
    node_coordinates_left, node_coordinates_right = node_coordinates
    solver_left, solver_right = solver
    compute_coefficients!(u_left, func, t, mesh.mesh_left, equations, solver_left, cache,
                          node_coordinates_left)
    compute_coefficients!(u_right, func, t, mesh.mesh_right, equations, solver_right, cache,
                          node_coordinates_right)
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

function calc_volume_integral!(du, u, mesh::OversetGridMesh, equations,
                               integral::Union{VolumeIntegralStrongForm,
                                               VolumeIntegralWeakForm},
                               solver, cache)
    u_left, u_right = u
    du_left, du_right = du
    f_all_left, f_all_right = cache.f_all
    volume_operator_left, volume_operator_right = cache.volume_operator
    solver_left, solver_right = solver
    calc_volume_integral!(du_left, u_left, mesh.mesh_left, equations,
                          integral, solver_left, volume_operator_left, f_all_left)
    calc_volume_integral!(du_right, u_right, mesh.mesh_right, equations,
                          integral, solver_right, volume_operator_right, f_all_right)
end

function create_cache_surface_flux_values(mesh::OversetGridMesh, equations, solver)
    solver_left, solver_right = solver
    surface_flux_values_left = create_cache_surface_flux_values(mesh.mesh_left, equations,
                                                                solver_left)
    surface_flux_values_right = create_cache_surface_flux_values(mesh.mesh_right, equations,
                                                                 solver_right)
    surface_flux_values = VectorOfArray([
                                            surface_flux_values_left,
                                            surface_flux_values_right
                                        ])
    return surface_flux_values
end

function calc_interface_flux!(surface_flux_values, u, mesh::OversetGridMesh, equations,
                              integral::AbstractSurfaceIntegral, solver, cache)
    u_left, u_right = u
    surface_flux_values_left, surface_flux_values_right = surface_flux_values
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    calc_interface_flux!(surface_flux_values_left, u_left, mesh_left,
                         equations, integral, solver_left, cache)
    calc_interface_flux!(surface_flux_values_right, u_right, mesh_right,
                         equations, integral, solver_right, cache)
end

function calc_boundary_flux!(surface_flux_values, u, t, boundary_conditions,
                             mesh::OversetGridMesh, equations,
                             integral::AbstractSurfaceIntegral, solver)
    u_left, u_right = u
    (; x_neg, x_pos) = boundary_conditions
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    surface_flux_values_left, surface_flux_values_right = surface_flux_values
    solver_left, solver_right = solver

    # Left boundary condition of left mesh
    u_ll = x_neg(u, xmin(mesh), t, mesh, equations, solver, true)
    u_rr = get_node_vars(u_left, equations, 1, 1)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
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
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
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
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values_right, f, equations, 1, 1)

    # Right boundary condition of right mesh
    u_ll = get_node_vars(u_right, equations, nnodes(solver_right, nelements(mesh_right)),
                         nelements(mesh_right))
    u_rr = x_pos(u, xmax(mesh), t, mesh, equations, solver, false)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values_right, f, equations, 2, nelements(mesh_right))
end

function calc_surface_integral!(du, u, mesh::OversetGridMesh, equations,
                                integral::SurfaceIntegralWeakForm,
                                solver, cache)
    du_left, du_right = du
    u_left, u_right = u
    surface_flux_values_left, surface_flux_values_right = cache.surface_flux_values
    solver_left, solver_right = solver
    surface_operator_left, surface_operator_right = cache.surface_operator
    calc_surface_integral!(du_left, u_left, mesh.mesh_left, equations,
                           integral, solver_left, surface_operator_left,
                           surface_flux_values_left)
    calc_surface_integral!(du_right, u_right, mesh.mesh_right, equations,
                           integral, solver_right, surface_operator_right,
                           surface_flux_values_right)
end

function calc_surface_integral!(du, u, mesh::OversetGridMesh, equations,
                                integral::SurfaceIntegralStrongForm,
                                solver, cache)
    du_left, du_right = du
    u_left, u_right = u
    surface_flux_values_left, surface_flux_values_right = cache.surface_flux_values
    solver_left, solver_right = solver
    surface_operator_left_left, surface_operator_right_left = cache.surface_operator_left
    surface_operator_left_right, surface_operator_right_right = cache.surface_operator_right
    calc_surface_integral!(du_left, u_left, mesh.mesh_left, equations,
                           integral, solver_left, surface_operator_left_left,
                           surface_operator_left_right, surface_flux_values_left)
    calc_surface_integral!(du_right, u_right, mesh.mesh_right, equations,
                           integral, solver_right, surface_operator_right_left,
                           surface_operator_right_right, surface_flux_values_right)
end

function apply_jacobian!(du, mesh::OversetGridMesh, equations, solver, cache)
    du_left, du_right = du
    jacobian_left, jacobian_right = cache.jacobian
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    apply_jacobian!(du_left, mesh_left, equations, solver_left, cache, jacobian_left)
    apply_jacobian!(du_right, mesh_right, equations, solver_right, cache, jacobian_right)
end

# This method is for integrating a vector quantity for all variables over the entire domain,
# such as the whole solution vector `u` (`VectorOfArray{T, 4, Vector{Array{T, 3}}}` for DG methods
# with same basis across elements and `VectorOfArray{T, 4, Vector{VectorOfArray{T, 3, Vector{Matrix{T}}}}}}`
# for `PerElementFDSBP`).
function PolynomialBases.integrate(func,
                                   u::Union{VectorOfArray{T, 4, Vector{Array{T, 3}}},
                                            VectorOfArray{T, 4,
                                                          Vector{VectorOfArray{T, 3,
                                                                               Vector{Matrix{T}}}}}},
                                   semi::SemidiscretizationOversetGrid) where {T}
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
                                   u::VectorOfArray{T, 3,
                                                    Vector{VectorOfArray{T, 2,
                                                                         Vector{Vector{T}}}}},
                                   semi::SemidiscretizationOversetGrid) where {T}
    u_left, u_right = u
    jacobian_left, jacobian_right = semi.cache.jacobian
    solver_left, solver_right = semi.solver
    l_left = left_overlap_element(semi.mesh)
    # TODO: This integrates the left domain only until the left boundary of the left overlap element,
    # i.e., we miss the integral from the left boundary of the left overlap element to b. Which part
    # should be integrated in the overlap region?
    integral_left = sum(integrate_on_element(func, u_left.u[element], solver_left, element,
                                             jacobian_left)
                        for element in 1:l_left)
    integral_right = sum(integrate_on_element(func, u_right.u[element], solver_right,
                                              element, jacobian_right)
                         for element in 1:nelements(semi.mesh.mesh_right))
    return integral_left + integral_right
end

function integrate_quantity!(quantity, func, u, semi::SemidiscretizationOversetGrid)
    mesh, equations, solver, _ = mesh_equations_solver_cache(semi)
    u_left, u_right = u
    quantity_left, quantity_right = quantity
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    compute_quantity!(quantity_left, func, u_left, mesh_left, equations, solver_left)
    compute_quantity!(quantity_right, func, u_right, mesh_right, equations, solver_right)
    integrate(quantity, semi)
end

function analyze(::typeof(entropy_timederivative), du, u, t,
                 semi::SemidiscretizationOversetGrid)
    mesh, equations, solver, _ = mesh_equations_solver_cache(semi)
    u_left, u_right = u
    du_left, du_right = du
    quantity = get_tmp_cache_scalar(semi)
    quantity_left, quantity_right = quantity
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    compute_quantity_timederivative!(quantity_left, cons2entropy, du_left, u_left,
                                     mesh_left, equations, solver_left)
    compute_quantity_timederivative!(quantity_right, cons2entropy, du_right, u_right,
                                     mesh_right, equations, solver_right)
    return integrate(quantity, semi)
end

function calc_error_norms(u, t, initial_condition, mesh::OversetGridMesh,
                          equations, solver, cache)
    u_left, u_right = u
    jacobian_left, jacobian_right = cache.jacobian
    node_coordinates_left, node_coordinates_right = cache.node_coordinates
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    l2_error_left, linf_error_left = calc_error_norms(u_left, t, initial_condition,
                                                      mesh_left, equations,
                                                      solver_left, cache, jacobian_left,
                                                      node_coordinates_left)
    l2_error_right, linf_error_right = calc_error_norms(u_right, t, initial_condition,
                                                        mesh_right, equations,
                                                        solver_right, cache, jacobian_right,
                                                        node_coordinates_right)
    return l2_error_left + l2_error_right, max(linf_error_left, linf_error_right)
end

function max_dt(u, t, mesh::OversetGridMesh, equations, solver, jacobian)
    u_left, u_right = u
    jacobian_left, jacobian_right = jacobian
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    solver_left, solver_right = solver
    max_scaled_speed_left = max_dt(u_left, t, mesh_left, equations, solver_left,
                                   jacobian_left)
    max_scaled_speed_right = max_dt(u_right, t, mesh_right, equations, solver_right,
                                    jacobian_right)
    return max(max_scaled_speed_left, max_scaled_speed_right)
end
