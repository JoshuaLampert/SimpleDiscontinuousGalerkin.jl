"""
    flux_central(u_ll, u_rr, orientation_or_normal_direction, equations::AbstractEquations)

The classical central numerical flux `f((u_ll) + f(u_rr)) / 2`. When this flux is
used as volume flux, the discretization is equivalent to the classical weak form
DG method (except floating point errors).
"""
@inline function flux_central(u_ll, u_rr, equations::AbstractEquations)
    # Calculate regular 1D fluxes
    f_ll = flux(u_ll, equations)
    f_rr = flux(u_rr, equations)

    # Average regular fluxes
    return 0.5f0 * (f_ll + f_rr)
end
