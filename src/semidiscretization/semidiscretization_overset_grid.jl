SemidiscretizationOversetGrid = Semidiscretization{<:OversetGridMesh}

# TODO: Allow different solvers for left and right mesh (e.g. for `PerElementBasis`)
function create_jacobian_and_node_coordinates(mesh::OversetGridMesh,
                                              solver::Union{DGSEM, FDSBP})
    jacobian_left, x_left = create_jacobian_and_node_coordinates(mesh.mesh_left, solver)
    jacobian_right, x_right = create_jacobian_and_node_coordinates(mesh.mesh_right, solver)
    return (jacobian_left, jacobian_right), (x_left, x_right)
end

function allocate_coefficients(mesh::OversetGridMesh, equations,
                               solver::Union{DGSEM, FDSBP})
    u_left = allocate_coefficients(mesh.mesh_left, equations, solver)
    u_right = allocate_coefficients(mesh.mesh_right, equations, solver)
    return VectorOfArray([u_left, u_right])
end

function compute_coefficients!(u, func, t, mesh::OversetGridMesh, equations, solver,
                               cache, node_coordinates)
    compute_coefficients!(u[1], func, t, mesh.mesh_left, equations, solver, cache,
                          node_coordinates[1])
    compute_coefficients!(u[2], func, t, mesh.mesh_right, equations, solver, cache,
                          node_coordinates[2])
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
