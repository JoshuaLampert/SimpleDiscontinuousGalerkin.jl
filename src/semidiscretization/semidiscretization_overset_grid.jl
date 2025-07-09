SemidiscretizationOversetGrid = Semidiscretization{<:OversetGridMesh}

# TODO: Allow different solvers for left and right mesh (e.g. for `PerElementBasis`)
function create_jacobian_and_node_coordinates(mesh::OversetGridMesh,
                                              solver::Union{DGSEM, FDSBP})
    jacobian_left, x_left = create_jacobian_and_node_coordinates(mesh.mesh_left, solver)
    jacobian_right, x_right = create_jacobian_and_node_coordinates(mesh.mesh_right, solver)
    return (jacobian_left, jacobian_right), (x_left, x_right)
end

function create_tmp_scalar(mesh::OversetGridMesh, solver)
    tmp_scalar_left = create_tmp_scalar(mesh.mesh_left, solver)
    tmp_scalar_right = create_tmp_scalar(mesh.mesh_right, solver)
    return VectorOfArray([tmp_scalar_left, tmp_scalar_right])
end

function allocate_coefficients(mesh::OversetGridMesh, equations,
                               solver::Union{DGSEM, FDSBP})
    u_left = allocate_coefficients(mesh.mesh_left, equations, solver)
    u_right = allocate_coefficients(mesh.mesh_right, equations, solver)
    return VectorOfArray([u_left, u_right])
end

function compute_coefficients!(u, func, t, mesh::OversetGridMesh, equations, solver,
                               cache, node_coordinates)
    u_left, u_right = u
    node_coordinates_left, node_coordinates_right = node_coordinates
    compute_coefficients!(u_left, func, t, mesh.mesh_left, equations, solver, cache,
                          node_coordinates_left)
    compute_coefficients!(u_right, func, t, mesh.mesh_right, equations, solver, cache,
                          node_coordinates_right)
end

function flat_grid(semi::SemidiscretizationOversetGrid)
    return vec(grid(semi)[1]), vec(grid(semi)[2])
end
# To avoid ambiguities
function flat_grid(::Semidiscretization{M, E, I, B, S}) where {M <: OversetGridMesh, E, I,
                                                               B, S <: PerElementFDSBP}
    return vec(grid(semi)[1]), vec(grid(semi)[2])
end
function get_variable(u, v, semi::SemidiscretizationOversetGrid)
    u_left, u_right = u
    return get_variable(u_left, v, semi.solver), get_variable(u_right, v, semi.solver)
end

function calc_volume_integral!(du, u, mesh::OversetGridMesh, equations,
                               integral::Union{VolumeIntegralStrongForm,
                                               VolumeIntegralWeakForm},
                               solver, cache)
    u_left, u_right = u
    du_left, du_right = du
    calc_volume_integral!(du_left, u_left, mesh.mesh_left, equations,
                          integral, solver, cache, cache.f_all[1])
    calc_volume_integral!(du_right, u_right, mesh.mesh_right, equations,
                          integral, solver, cache, cache.f_all[2])
end

function create_cache_surface_flux_values(mesh::OversetGridMesh, equations, solver)
    surface_flux_values_left = create_cache_surface_flux_values(mesh.mesh_left, equations,
                                                                solver)
    surface_flux_values_right = create_cache_surface_flux_values(mesh.mesh_right, equations,
                                                                 solver)
    surface_flux_values = VectorOfArray([
                                            surface_flux_values_left,
                                            surface_flux_values_right
                                        ])
    return surface_flux_values
end

function calc_interface_flux!(surface_flux_values, u, mesh::OversetGridMesh, equations,
                              integral::AbstractSurfaceIntegral, solver, cache)
    u_left, u_right = u
    calc_interface_flux!(surface_flux_values[1], u_left, mesh.mesh_left,
                         equations, integral, solver, cache)
    calc_interface_flux!(surface_flux_values[2], u_right, mesh.mesh_right,
                         equations, integral, solver, cache)
end

function calc_boundary_flux!(surface_flux_values, u, t, boundary_conditions,
                             mesh::OversetGridMesh, equations,
                             integral::AbstractSurfaceIntegral, solver)
    u_left, u_right = u
    (; x_neg, x_pos) = boundary_conditions
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right

    # Left boundary condition of left mesh
    u_ll = x_neg(u, xmin(mesh), t, mesh, equations, solver, true)
    u_rr = get_node_vars(u_left, equations, 1, 1)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values[1], f, equations, 1, 1)

    # Right boundary condition of left mesh
    u_ll = get_node_vars(u_left, equations, nnodes(solver, nelements(mesh_left)),
                         nelements(mesh_left))
    # TODO: Make this more efficient by computing e_M_right and e_M_left only once
    # during initialization and storing it in the cache
    l_right = right_overlap_element(mesh)
    c = xmax(mesh_left)
    D = get_basis(solver, l_right)
    linear_map(x, a, b, c, d) = c + (x - a) / (b - a) * (d - c)
    xlr_L = left_element_boundary(mesh_right, l_right)
    xlr_R = left_element_boundary(mesh_right, l_right + 1)
    c_mapped = linear_map(c, xlr_L, xlr_R, first(grid(D)), last(grid(D)))
    u_rr = zeros(real(solver), nvariables(equations))
    e_M_right = interpolation_operator(c_mapped, D)
    for v in eachvariable(equations)
        values = u_right[v, :, l_right]
        u_rr[v] = e_M_right' * values
    end
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values[1], f, equations, 2, nelements(mesh_left))

    # Left boundary condition of right mesh
    l_left = left_overlap_element(mesh)
    b = xmin(mesh_right)
    D = get_basis(solver, l_left)
    xll_L = left_element_boundary(mesh_left, l_left)
    xll_R = left_element_boundary(mesh_left, l_left + 1)
    b_mapped = linear_map(b, xll_L, xll_R, first(grid(D)), last(grid(D)))
    u_ll = zeros(real(solver), nvariables(equations))
    e_M_left = interpolation_operator(b_mapped, D)
    for v in eachvariable(equations)
        values = u_left[v, :, l_left]
        u_ll[v] = e_M_left' * values
    end
    u_rr = get_node_vars(u_right, equations, 1, 1)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values[2], f, equations, 1, 1)

    # Right boundary condition of right mesh
    u_ll = get_node_vars(u_right, equations, nnodes(solver, nelements(mesh_right)),
                         nelements(mesh_right))
    u_rr = x_pos(u, xmax(mesh), t, mesh, equations, solver, false)
    f = integral.surface_flux_boundary(u_ll, u_rr, equations)
    set_node_vars!(surface_flux_values[2], f, equations, 2, nelements(mesh_right))
end

function calc_surface_integral!(du, u, mesh::OversetGridMesh, equations,
                                integral::Union{SurfaceIntegralStrongForm,
                                                SurfaceIntegralWeakForm},
                                solver, cache)
    du_left, du_right = du
    u_left, u_right = u
    calc_surface_integral!(du_left, u_left, mesh.mesh_left, equations,
                           integral, solver, cache, cache.surface_flux_values[1])
    calc_surface_integral!(du_right, u_right, mesh.mesh_right, equations,
                           integral, solver, cache, cache.surface_flux_values[2])
end

function apply_jacobian!(du, mesh::OversetGridMesh, equations, solver, cache)
    du_left, du_right = du
    apply_jacobian!(du_left, mesh.mesh_left, equations, solver, cache, cache.jacobian[1])
    apply_jacobian!(du_right, mesh.mesh_right, equations, solver, cache, cache.jacobian[2])
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
    display(typeof(u))
    u_left, u_right = u
    jacobian_left, jacobian_right = semi.cache.jacobian
    l_left = left_overlap_element(semi.mesh)
    integral_left = sum(integrate_on_element(func, u_left.u[element], semi, element,
                                             jacobian_left)
                        for element in 1:l_left)
    integral_right = sum(integrate_on_element(func, u_right.u[element], semi, element,
                                              jacobian_right)
                         for element in 1:nelements(semi.mesh.mesh_right))
    return integral_left + integral_right
end

function integrate_quantity!(quantity, func, u, semi::SemidiscretizationOversetGrid)
    mesh, equations, solver, _ = mesh_equations_solver_cache(semi)
    u_left, u_right = u
    quantity_left, quantity_right = quantity
    mesh_left, mesh_right = mesh.mesh_left, mesh.mesh_right
    compute_quantity!(quantity_left, func, u_left, mesh_left, equations, solver)
    compute_quantity!(quantity_right, func, u_right, mesh_right, equations, solver)
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
    compute_quantity_timederivative!(quantity_left, cons2entropy, du_left, u_left,
                                     mesh_left, equations, solver)
    compute_quantity_timederivative!(quantity_right, cons2entropy, du_right, u_right,
                                     mesh_right, equations, solver)
    return integrate(quantity, semi)
end

function calc_error_norms(u, t, initial_condition, mesh::OversetGridMesh,
                          equations, solver, cache)
    u_left, u_right = u
    jacobian_left, jacobian_right = cache.jacobian
    node_coordinates_left, node_coordinates_right = cache.node_coordinates
    l2_error_left, linf_error_left = calc_error_norms(u_left, t, initial_condition,
                                                      mesh.mesh_left, equations,
                                                      solver, cache, jacobian_left,
                                                      node_coordinates_left)
    l2_error_right, linf_error_right = calc_error_norms(u_right, t, initial_condition,
                                                        mesh.mesh_right, equations,
                                                        solver, cache, jacobian_right,
                                                        node_coordinates_right)
    return l2_error_left + l2_error_right, max(linf_error_left, linf_error_right)
end
