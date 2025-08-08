@doc raw"""
    CompressibleEulerEquations1D(gamma)

The compressible Euler equations
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho e
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
\rho v_1 \\ \rho v_1^2 + p \\ (\rho e +p) v_1
\end{pmatrix}
=
\begin{pmatrix}
0 \\ 0 \\ 0
\end{pmatrix}
```
for an ideal gas with ratio of specific heats `gamma` in one space dimension.
Here, ``\rho`` is the density, ``v_1`` the velocity, ``e`` the specific total energy **rather than** specific internal energy, and
```math
p = (\gamma - 1) \left( \rho e - \frac{1}{2} \rho v_1^2 \right)
```
the pressure.
"""
struct CompressibleEulerEquations1D{RealT <: Real} <: AbstractEquations{1, 3}
    gamma::RealT               # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications

    function CompressibleEulerEquations1D(gamma)
        γ, inv_gamma_minus_one = promote(gamma, inv(gamma - 1))
        new{typeof(γ)}(γ, inv_gamma_minus_one)
    end
end

function varnames(::typeof(cons2cons), ::CompressibleEulerEquations1D)
    ("rho", "rho_v1", "rho_e")
end
varnames(::typeof(cons2prim), ::CompressibleEulerEquations1D) = ("rho", "v1", "p")

"""
    initial_condition_density_wave(x, t, equations::CompressibleEulerEquations1D)

A sine wave in the density with constant velocity and pressure; reduces the
compressible Euler equations to the linear advection equations.
This setup is the test case for stability of EC fluxes from the paper
- Gregor J. Gassner, Magnus Svärd, Florian J. Hindenlang (2020)
  Stability issues of entropy-stable and/or split-form high-order schemes
  [arXiv: 2007.09026](https://arxiv.org/abs/2007.09026)
with the following parameters
- domain [-1, 1]
- mesh = 4x4
- polydeg = 5
"""
function initial_condition_density_wave(x, t, equations::CompressibleEulerEquations1D)
    RealT = eltype(x)
    v1 = convert(RealT, 0.1)
    rho = 1 + convert(RealT, 0.98) * sinpi(2 * (x - t * v1))
    rho_v1 = rho * v1
    p = 20
    rho_e = p / (equations.gamma - 1) + 0.5f0 * rho * v1^2
    return SVector(rho, rho_v1, rho_e)
end

"""
    initial_condition_convergence_test(x, t, equations::CompressibleEulerEquations1D)

A smooth initial condition used for convergence tests in combination with
[`source_terms_convergence_test`](@ref)
(and [`BoundaryConditionDirichlet(initial_condition_convergence_test)`](@ref) in non-periodic domains).
"""
function initial_condition_convergence_test(x, t,
                                            equations::CompressibleEulerEquations1D)
    RealT = eltype(x)
    c = 2
    A = convert(RealT, 0.1)
    L = 2
    f = 1.0f0 / L
    ω = 2 * convert(RealT, pi) * f
    ini = c + A * sin(ω * (x - t))

    rho = ini
    rho_v1 = ini
    rho_e = ini^2

    return SVector(rho, rho_v1, rho_e)
end

"""
    source_terms_convergence_test(u, x, t, equations::CompressibleEulerEquations1D)

Source terms used for convergence tests in combination with
[`initial_condition_convergence_test`](@ref)
(and [`BoundaryConditionDirichlet(initial_condition_convergence_test)`](@ref) in non-periodic domains).
"""
@inline function source_terms_convergence_test(u, x, t,
                                               equations::CompressibleEulerEquations1D)
    # Same settings as in `initial_condition`
    RealT = eltype(u)
    c = 2
    A = convert(RealT, 0.1)
    L = 2
    f = 1.0f0 / L
    ω = 2 * convert(RealT, pi) * f
    γ = equations.gamma

    si, co = sincos(ω * (x - t))
    rho = c + A * si
    rho_x = ω * A * co

    # Note that d/dt rho = -d/dx rho.
    # This yields du2 = du3 = d/dx p (derivative of pressure).
    # Other terms vanish because of v = 1.
    du1 = 0
    du2 = rho_x * (2 * rho - 0.5f0) * (γ - 1)
    du3 = du2

    return SVector(du1, du2, du3)
end

"""
    initial_condition_weak_blast_wave(x, t, equations::CompressibleEulerEquations1D)

A weak blast wave taken from
- Sebastian Hennemann, Gregor J. Gassner (2020)
  A provably entropy stable subcell shock capturing approach for high order split form DG
  [arXiv: 2008.12044](https://arxiv.org/abs/2008.12044)
"""
function initial_condition_weak_blast_wave(x, t,
                                           equations::CompressibleEulerEquations1D)
    # From Hennemann & Gassner JCP paper 2020 (Sec. 6.3)
    # Set up polar coordinates
    RealT = eltype(x)
    # We often have an interface at x = 0, so we shift the center
    # to avoid having a discontinuity at an interface.
    inicenter = 0.001
    x_norm = x - inicenter
    r = abs(x_norm)
    # The following code is equivalent to
    # phi = atan(0.0, x_norm)
    # cos_phi = cos(phi)
    # in 1D but faster
    cos_phi = x_norm > 0 ? 1 : -1

    # Calculate primitive variables
    rho = r > 0.5f0 ? one(RealT) : convert(RealT, 1.1691)
    v1 = r > 0.5f0 ? zero(RealT) : convert(RealT, 0.1882) * cos_phi
    p = r > 0.5f0 ? one(RealT) : convert(RealT, 1.245)

    return prim2cons(SVector(rho, v1, p), equations)
end

# Calculate 1D flux for a single point
@inline function flux(u, equations::CompressibleEulerEquations1D)
    rho, rho_v1, rho_e = u
    v1 = rho_v1 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho_v1 * v1)
    # Ignore orientation since it is always "1" in 1D
    f1 = rho_v1
    f2 = rho_v1 * v1 + p
    f3 = (rho_e + p) * v1
    return SVector(f1, f2, f3)
end

"""
    flux_kennedy_gruber(u_ll, u_rr, equations::CompressibleEulerEquations1D)

Kinetic energy preserving two-point flux by
- Kennedy and Gruber (2008)
  Reduced aliasing formulations of the convective terms within the
  Navier-Stokes equations for a compressible fluid
  [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function flux_kennedy_gruber(u_ll, u_rr, equations::CompressibleEulerEquations1D)
    # Unpack left and right state
    rho_e_ll = last(u_ll)
    rho_e_rr = last(u_rr)
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Ignore orientation since it is always "1" in 1D
    f1 = rho_avg * v1_avg
    f2 = rho_avg * v1_avg * v1_avg + p_avg
    f3 = (rho_avg * e_avg + p_avg) * v1_avg

    return SVector(f1, f2, f3)
end

"""
    flux_ranocha(u_ll, u_rr, equations::CompressibleEulerEquations1D)

Entropy conserving and kinetic energy preserving two-point flux by
- Hendrik Ranocha (2018)
  Generalised Summation-by-Parts Operators and Entropy Stability of Numerical Methods
  for Hyperbolic Balance Laws
  [PhD thesis, TU Braunschweig](https://cuvillier.de/en/shop/publications/7743)
See also
- Hendrik Ranocha (2020)
  Entropy Conserving and Kinetic Energy Preserving Numerical Methods for
  the Euler Equations Using Summation-by-Parts Operators
  [Proceedings of ICOSAHOM 2018](https://doi.org/10.1007/978-3-030-39647-3_42)
"""
@inline function flux_ranocha(u_ll, u_rr, equations::CompressibleEulerEquations1D)
    # Unpack left and right state
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
    # in exact arithmetic since
    #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
    #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    velocity_square_avg = 0.5f0 * (v1_ll * v1_rr)

    # Calculate fluxes
    # Ignore orientation since it is always "1" in 1D
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = f1 * (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
         0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)

    return SVector(f1, f2, f3)
end

# Based on https://gist.github.com/ketch/08ce0845da0c8f3fa9ff
@inline function flux_godunov(u_ll, u_rr, equations::CompressibleEulerEquations1D)
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    gamma = equations.gamma

    c_ll = sqrt(gamma * p_ll / rho_ll)
    c_rr = sqrt(gamma * p_rr / rho_rr)

    alpha = (gamma - 1) / (2 * gamma)
    beta = (gamma + 1) / (gamma - 1)

    function phi_l(p)
        # Hugoniot locus 1
        if p >= p_ll
            v1_ll +
            2 * c_ll / sqrt(2 * gamma * (gamma - 1.0)) *
            ((1 - p / p_ll) / sqrt(1 + beta * p / p_ll))
        else # Integral curve 1
            v1_ll +
            2 * c_ll / (gamma - 1.0) * (1 - (p / p_ll)^alpha)
        end
    end
    function phi_r(p)
        # Hugoniot locus 3
        if p >= p_rr
            v1_rr -
            2 * c_rr / sqrt(2 * gamma * (gamma - 1.0)) *
            ((1 - p / p_rr) / sqrt(1 + beta * p / p_rr))
        else # Integral curve 3
            v1_rr -
            2 * c_rr / (gamma - 1.0) * (1 - (p / p_rr)^alpha)
        end
    end

    phi(p) = phi_l(p) - phi_r(p)
    p_star = find_zero(phi, (0.5f0 * min(p_ll, p_rr), 1.5f0 * max(p_ll, p_rr)),
                       AlefeldPotraShi())
    v1_star = phi_l(p_star)

    rho_ll_star = (p_star / p_ll)^(1.0 / gamma) * rho_ll
    rho_rr_star = (p_star / p_rr)^(1.0 / gamma) * rho_rr

    w3 = v1_star
    if p_star > p_ll
        w1 = (rho_ll * v1_ll - rho_ll_star * v1_star) / (rho_ll - rho_ll_star)
        w2 = w1
    else
        c_ll_star = sqrt(gamma * p_star / rho_ll_star)
        w1 = v1_ll - c_ll
        w2 = v1_star - c_ll_star
    end
    if p_star > p_rr
        w5 = (rho_rr * v1_rr - rho_rr_star * v1_star) / (rho_rr - rho_rr_star)
        w4 = w5
    else
        c_rr_star = sqrt(gamma * p_star / rho_rr_star)
        w4 = v1_star + c_rr_star
        w5 = v1_rr + c_rr
    end

    xi = 0.0 # x/t, but Godunov flux uses the exact Riemann solution at x/t = 0
    v1_1 = ((gamma - 1.0) * v1_ll + 2 * (c_ll + xi)) / (gamma + 1.0)
    v1_3 = ((gamma - 1.0) * v1_rr - 2 * (c_rr - xi)) / (gamma + 1.0)
    rho_1 = (rho_ll^gamma * (v1_1 - xi)^2 / (gamma * p_ll))^equations.inv_gamma_minus_one
    rho_3 = (rho_rr^gamma * (xi - v1_3)^2 / (gamma * p_rr))^equations.inv_gamma_minus_one
    p_1 = p_ll * (rho_1 / rho_ll)^gamma
    p_3 = p_rr * (rho_3 / rho_rr)^gamma

    if xi <= w1
        prim_out = SVector(rho_ll, v1_ll, p_ll)
    elseif w1 < xi <= w2
        prim_out = SVector(rho_1, v1_1, p_1)
    elseif w2 < xi <= w3
        prim_out = SVector(rho_ll_star, v1_star, p_star)
    elseif w3 < xi <= w4
        prim_out = SVector(rho_rr_star, v1_star, p_star)
    elseif w4 < xi <= w5
        prim_out = SVector(rho_3, v1_3, p_3)
    else # xi > w5
        prim_out = SVector(rho_rr, v1_rr, p_rr)
    end
    return flux(prim2cons(prim_out, equations), equations)
end

@inline function max_abs_speed(u_ll, u_rr, equations::CompressibleEulerEquations1D)
    rho_ll, rho_v1_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_e_rr = u_rr

    # Calculate primitive variables and speed of sound
    v1_ll = rho_v1_ll / rho_ll
    v_mag_ll = abs(v1_ll)
    p_ll = (equations.gamma - 1) * (rho_e_ll - 0.5f0 * rho_ll * v_mag_ll^2)
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    v1_rr = rho_v1_rr / rho_rr
    v_mag_rr = abs(v1_rr)
    p_rr = (equations.gamma - 1) * (rho_e_rr - 0.5f0 * rho_rr * v_mag_rr^2)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(v_mag_ll + c_ll, v_mag_rr + c_rr)
end

# For the HLL flux
@inline function min_max_speed_davis(u_ll, u_rr, equations::CompressibleEulerEquations1D)
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    λ_min = min(v1_ll - c_ll, v1_rr - c_rr)
    λ_max = max(v1_ll + c_ll, v1_rr + c_rr)

    return λ_min, λ_max
end

# Convert conservative variables to primitive
@inline function cons2prim(u, equations::CompressibleEulerEquations1D)
    rho, rho_v1, rho_e = u

    v1 = rho_v1 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho_v1 * v1)

    return SVector(rho, v1, p)
end

# Convert conservative variables to entropy
@inline function cons2entropy(u, equations::CompressibleEulerEquations1D)
    rho, rho_v1, rho_e = u

    v1 = rho_v1 / rho
    v_square = v1^2
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * v_square)
    s = log(p) - equations.gamma * log(rho)
    rho_p = rho / p

    w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
         0.5f0 * rho_p * v_square
    w2 = rho_p * v1
    w3 = -rho_p

    return SVector(w1, w2, w3)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::CompressibleEulerEquations1D)
    rho, v1, p = prim
    rho_v1 = rho * v1
    rho_e = p * equations.inv_gamma_minus_one + 0.5f0 * (rho_v1 * v1)
    return SVector(rho, rho_v1, rho_e)
end

@inline function density(u, equations::CompressibleEulerEquations1D)
    rho = u[1]
    return rho
end

@inline function velocity(u, equations::CompressibleEulerEquations1D)
    rho = u[1]
    v1 = u[2] / rho
    return v1
end

@inline function pressure(u, equations::CompressibleEulerEquations1D)
    rho, rho_v1, rho_e = u
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1^2) / rho)
    return p
end

@inline function density_pressure(u, equations::CompressibleEulerEquations1D)
    rho, rho_v1, rho_e = u
    rho_times_p = (equations.gamma - 1) * (rho * rho_e - 0.5f0 * (rho_v1^2))
    return rho_times_p
end

# Calculate thermodynamic entropy for a conservative state `cons`
@inline function entropy_thermodynamic(cons, equations::CompressibleEulerEquations1D)
    # Pressure
    p = (equations.gamma - 1) * (cons[3] - 0.5f0 * (cons[2]^2) / cons[1])

    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])

    return s
end

# Calculate mathematical entropy for a conservative state `cons`
@inline function entropy_math(cons, equations::CompressibleEulerEquations1D)
    # Mathematical entropy
    S = -entropy_thermodynamic(cons, equations) * cons[1] *
        equations.inv_gamma_minus_one

    return S
end

# Default entropy is the mathematical entropy
@inline function entropy(cons, equations::CompressibleEulerEquations1D)
    entropy_math(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline energy_total(cons, ::CompressibleEulerEquations1D) = cons[3]
