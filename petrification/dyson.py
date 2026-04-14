"""
Dyson equation as a fixed-point problem, with alpha-transform relaxation.

The Dyson equation for the dressed propagator is:
    G = G₀ + G₀ Σ[G] G
or equivalently:
    G(ω) = 1 / (ω - ε₀ - Σ[G](ω))

This is literally a fixed-point equation G = F(G). Standard perturbation
theory iterates F naively; for strong coupling this diverges. Alpha-
relaxation replaces the iteration with:
    G_{n+1} = α·F(G_n) + (1-α)·G_n

This module provides scalar, matrix, and lattice Dyson equation solvers
with alpha-relaxation, including frequency-dependent α(ω).
"""

import numpy as np


# ============================================================
# Scalar Dyson equation: G = 1/(ω - ε₀ - Σ(G))
# ============================================================

def scalar_dyson_iterate(omega, eps0, sigma_func, G0, alpha=1.0,
                         n_iter=200, tol=1e-12):
    """
    Iterate the scalar Dyson equation G = 1/(ω - ε₀ - Σ(G)) with
    alpha-relaxation.

    Parameters
    ----------
    omega : float or complex
        Frequency.
    eps0 : float
        Bare energy level.
    sigma_func : callable
        Self-energy functional Σ(G). Takes a complex number, returns complex.
    G0 : complex
        Initial guess for G.
    alpha : float or callable
        Relaxation parameter. If callable, alpha(omega, G, iteration).
    n_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance on |G_{n+1} - G_n|.

    Returns
    -------
    G_history : ndarray of complex
        G at each iteration.
    converged : bool
        Whether iteration converged within tol.
    """
    G = complex(G0)
    history = [G]

    for i in range(n_iter):
        sigma = sigma_func(G)
        denom = omega - eps0 - sigma
        if abs(denom) < 1e-30:
            # Dyson map has a singularity
            history.append(G)
            break

        F_G = 1.0 / denom

        a = alpha(omega, G, i) if callable(alpha) else alpha
        G_new = a * F_G + (1 - a) * G

        history.append(G_new)

        if abs(G_new - G) < tol:
            return np.array(history), True

        G = G_new

    return np.array(history), False


def scalar_dyson_exact(omega, eps0, U):
    """
    Exact solution of G = 1/(ω - ε₀ - U²G).

    This is a quadratic: U²G² - (ω - ε₀)G + 1 = 0
    G = [(ω - ε₀) ± √((ω - ε₀)² - 4U²)] / (2U²)

    We pick the physical root (the one that reduces to G₀ = 1/(ω - ε₀)
    as U → 0, i.e., the root with the correct sign of Im(G) for
    retarded propagator).

    Parameters
    ----------
    omega : float or complex
        Frequency (add +iη for retarded).
    eps0 : float
        Bare level.
    U : float
        Coupling strength.

    Returns
    -------
    G : complex
        Exact dressed propagator.
    """
    if abs(U) < 1e-15:
        return 1.0 / (omega - eps0)

    disc = (omega - eps0)**2 - 4 * U**2
    sqrt_disc = np.sqrt(complex(disc))

    G_plus = ((omega - eps0) + sqrt_disc) / (2 * U**2)
    G_minus = ((omega - eps0) - sqrt_disc) / (2 * U**2)

    # Physical root: for retarded propagator (ω + iη), Im(G) < 0
    if np.imag(omega) > 0:
        # Pick root with negative imaginary part
        return G_minus if np.imag(G_minus) < np.imag(G_plus) else G_plus
    else:
        # For real ω, pick root continuous with G₀ as U → 0
        G0 = 1.0 / (omega - eps0) if abs(omega - eps0) > 1e-15 else 1e15
        if abs(G_plus - G0) < abs(G_minus - G0):
            return G_plus
        return G_minus


def scalar_dyson_stability(omega, eps0, U):
    """
    Compute |F'(G*)| for the scalar Dyson map F(G) = 1/(ω - ε₀ - U²G).

    F'(G) = U² / (ω - ε₀ - U²G)² = U² G²

    The fixed point is stable under naive iteration iff |F'(G*)| < 1.

    Returns
    -------
    stability : float
        |F'(G*)| — the Lyapunov multiplier.
    G_star : complex
        The fixed point.
    alpha_opt : float
        Optimal alpha = 1/(1 - F'(G*)) for zero derivative (when real).
    """
    G_star = scalar_dyson_exact(omega, eps0, U)
    F_prime = U**2 * G_star**2
    stability = abs(F_prime)

    # Optimal alpha (from our alpha-transform theory)
    if abs(1 - F_prime) > 1e-10:
        alpha_opt = 1.0 / (1.0 - F_prime)
    else:
        alpha_opt = complex(float('inf'))

    return stability, G_star, alpha_opt


# ============================================================
# Spectral function Dyson equation (frequency-resolved)
# ============================================================

def spectral_dyson_scan(omega_grid, eps0, sigma_func, alpha=1.0,
                        n_iter=200, tol=1e-10):
    """
    Solve the Dyson equation at each frequency independently.

    Parameters
    ----------
    omega_grid : ndarray of complex
        Frequencies (typically ω + iη).
    eps0 : float
        Bare level.
    sigma_func : callable
        Σ(G) — self-energy as function of G.
    alpha : float or callable
        If callable, alpha(omega, G, iteration) — can be frequency-dependent.

    Returns
    -------
    G_converged : ndarray of complex
        Converged G(ω) at each frequency.
    n_converged : int
        Number of frequencies that converged.
    iterations : ndarray of int
        Iteration count at each frequency.
    """
    n_omega = len(omega_grid)
    G_out = np.zeros(n_omega, dtype=complex)
    converged_count = 0
    iter_counts = np.zeros(n_omega, dtype=int)

    for j, omega in enumerate(omega_grid):
        # Initial guess: bare propagator
        G0 = 1.0 / (omega - eps0)
        history, conv = scalar_dyson_iterate(omega, eps0, sigma_func, G0,
                                              alpha=alpha, n_iter=n_iter,
                                              tol=tol)
        G_out[j] = history[-1]
        iter_counts[j] = len(history) - 1
        if conv:
            converged_count += 1

    return G_out, converged_count, iter_counts


# ============================================================
# Lattice Dyson equation: G(ω) = [ω I - H₀ - Σ(G)]⁻¹
# ============================================================

def lattice_dyson_iterate(omega, H0, sigma_func, G0, alpha=1.0,
                          n_iter=200, tol=1e-10):
    """
    Iterate the matrix Dyson equation G = (ωI - H₀ - Σ[G])⁻¹
    with alpha-relaxation.

    Parameters
    ----------
    omega : complex
        Frequency.
    H0 : ndarray, shape (n, n)
        Bare Hamiltonian.
    sigma_func : callable
        Matrix self-energy Σ(G) → ndarray (n, n).
    G0 : ndarray, shape (n, n)
        Initial guess.
    alpha : float
        Relaxation parameter.
    n_iter : int
    tol : float

    Returns
    -------
    G_final : ndarray
        Converged Green's function matrix.
    residuals : list of float
        ||G_{n+1} - G_n|| at each iteration.
    converged : bool
    """
    n = H0.shape[0]
    G = G0.copy().astype(complex)
    residuals = []

    for i in range(n_iter):
        sigma = sigma_func(G)
        A = omega * np.eye(n) - H0 - sigma
        try:
            F_G = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            residuals.append(float('inf'))
            break

        G_new = alpha * F_G + (1 - alpha) * G
        res = np.linalg.norm(G_new - G)
        residuals.append(res)

        if res < tol:
            return G_new, residuals, True

        G = G_new

    return G, residuals, False


# ============================================================
# Self-energy models
# ============================================================

def quadratic_self_energy(U):
    """Σ(G) = U²·G — simplest self-consistent model."""
    return lambda G: U**2 * G


def cubic_self_energy(U, V):
    """Σ(G) = U²·G + V·G² — nonlinear self-energy."""
    return lambda G: U**2 * G + V * G**2


def hubbard_atom_self_energy(U, beta=10.0):
    """
    Self-energy for the single-site Hubbard model (atomic limit).

    In the atomic limit, the exact self-energy for the half-filled
    Hubbard atom is Σ(ω) = U²/(4ω), independent of G. But for
    testing purposes, we use a G-dependent approximation:
    Σ(G) ≈ U²·G (Hartree-Fock-like, second-order in U).
    """
    return lambda G: U**2 * G
