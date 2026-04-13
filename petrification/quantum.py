"""
Quantum mechanics tools: Hamiltonian discretization, Riccati-Bloch integration,
spectral scanning via shooting method, and alpha-relaxed Riccati iteration.

Bridges the gap between fixed-point dynamics and quantum eigenvalue problems
through Turbiner's nonlinearization procedure.
"""

import numpy as np
from scipy.linalg import eigh, solve
from scipy.integrate import solve_ivp


# ============================================================
# Hamiltonian construction
# ============================================================

def discretize_hamiltonian(V_func, x_grid, hbar=1.0, m=1.0):
    """
    Build the Hamiltonian matrix H = -(hbar^2/2m) d^2/dx^2 + V(x)
    using second-order finite differences.

    Parameters
    ----------
    V_func : callable
        Potential V(x).
    x_grid : ndarray
        Spatial grid points.
    hbar, m : float
        Physical constants (default: natural units).

    Returns
    -------
    H : ndarray
        (N x N) Hamiltonian matrix.
    """
    n = len(x_grid)
    h = x_grid[1] - x_grid[0]
    coeff = hbar**2 / (2.0 * m * h**2)

    # Kinetic energy via finite differences
    diag = np.full(n, 2.0 * coeff)
    off = np.full(n - 1, -coeff)
    T = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)

    # Potential energy
    V = np.diag(V_func(x_grid))

    return T + V


def solve_eigenstates(V_func, x_grid, n_states=None, **kwargs):
    """
    Solve for eigenvalues and eigenvectors of H = T + V.

    Returns
    -------
    eigenvalues : ndarray
    eigenvectors : ndarray
        Column k is the k-th eigenstate.
    """
    H = discretize_hamiltonian(V_func, x_grid, **kwargs)
    eigenvalues, eigenvectors = eigh(H)
    if n_states is not None:
        eigenvalues = eigenvalues[:n_states]
        eigenvectors = eigenvectors[:, :n_states]
    return eigenvalues, eigenvectors


# ============================================================
# Riccati-Bloch nonlinearization (Turbiner's method)
# ============================================================

def riccati_solve(V_func, E, x_span, y0=0.0, max_y=100.0,
                  max_step=0.01, rtol=1e-8, atol=1e-10):
    """
    Integrate the Riccati-Bloch equation

        y'(x) = y(x)^2 - 2*V(x) + 2*E

    where y = -psi'/psi is the logarithmic derivative.

    From -hbar^2/(2m) psi'' + V psi = E psi (with hbar=m=1):
        psi = exp(-int y dx)
        y' = y^2 - 2V + 2E

    At eigenvalues, y(x) remains bounded across the domain.
    Away from eigenvalues, y(x) diverges (Riccati singularity).

    Parameters
    ----------
    V_func : callable
        Potential function V(x).
    E : float
        Trial energy.
    x_span : tuple
        Integration interval (x_start, x_end).
    y0 : float
        Initial value of y.
    max_y : float
        Divergence cutoff.

    Returns
    -------
    sol : OdeSolution
        Solution object from scipy.integrate.solve_ivp.
    """
    def rhs(x, y):
        return [y[0]**2 - 2.0 * V_func(x) + 2.0 * E]

    def divergence_event(x, y):
        return max_y - abs(y[0])
    divergence_event.terminal = True
    divergence_event.direction = -1

    sol = solve_ivp(rhs, x_span, [y0], method='RK45',
                    events=[divergence_event],
                    max_step=max_step, dense_output=True,
                    rtol=rtol, atol=atol)
    return sol


def riccati_alpha(V_func, E, x_grid, y0=0.0, alpha=1.0, max_y=1e4):
    """
    Alpha-relaxed Riccati-Bloch iteration (forward Euler).

        y_{n+1} = alpha * [y_n + h*(y_n^2 - 2V + 2E)] + (1-alpha) * y_n

    This is the discrete analog of your alpha-transform applied to
    Turbiner's nonlinearization. Alpha controls the stability basin:
    different alpha values can stabilize trajectories for excited states.

    Returns
    -------
    trajectory : ndarray
        y(x) values on the grid.
    """
    h = x_grid[1] - x_grid[0]
    y = float(y0)
    trajectory = np.zeros(len(x_grid))
    trajectory[0] = y

    for i in range(1, len(x_grid)):
        x = x_grid[i - 1]
        y_rhs = y**2 - 2.0 * V_func(x) + 2.0 * E
        y_euler = y + h * y_rhs
        y_new = alpha * y_euler + (1 - alpha) * y

        if abs(y_new) > max_y:
            trajectory[i:] = np.sign(y_new) * max_y
            break
        y = y_new
        trajectory[i] = y

    return trajectory


# ============================================================
# Spectral scanning (shooting method)
# ============================================================

def spectral_scan(V_func, E_range, x_span=(0.0, 7.0), parity="even"):
    """
    For each trial energy E, integrate the Schrödinger equation
    as a shooting problem and record psi(x_max).

    At eigenvalues: psi(x_max) -> 0 (boundary condition satisfied).
    Away from eigenvalues: psi(x_max) diverges.

    This casting is itself a fixed-point problem: find E such that
    the shooting map S(E) = psi(x_max; E) satisfies S(E) = 0.

    Parameters
    ----------
    V_func : callable
        Potential function.
    E_range : ndarray
        Trial energies.
    x_span : tuple
        Integration domain.
    parity : str
        "even" (psi(0)=1, psi'(0)=0) or "odd" (psi(0)=0, psi'(0)=1).

    Returns
    -------
    psi_end : ndarray
        psi(x_max) for each trial energy.
    """
    if parity == "even":
        ic = [1.0, 0.0]
    elif parity == "odd":
        ic = [0.0, 1.0]
    else:
        raise ValueError(f"parity must be 'even' or 'odd', got '{parity}'")

    results = np.zeros(len(E_range))

    for j, E in enumerate(E_range):
        def schrodinger(x, state):
            psi, dpsi = state
            return [dpsi, 2.0 * (V_func(x) - E) * psi]

        sol = solve_ivp(schrodinger, list(x_span), ic,
                        method='RK45', max_step=0.02,
                        rtol=1e-10, atol=1e-12)
        results[j] = sol.y[0, -1]

    return results


def detect_eigenvalues(V_func, E_range, x_span=(0.0, 7.0), order=20):
    """
    Detect eigenvalues by finding minima of |psi(x_max)| for both parities.

    Returns
    -------
    E_detected : ndarray
        Detected eigenvalue positions.
    """
    from scipy.signal import argrelmin

    psi_even = spectral_scan(V_func, E_range, x_span, parity="even")
    psi_odd = spectral_scan(V_func, E_range, x_span, parity="odd")

    peaks_even = argrelmin(np.abs(psi_even), order=order)[0]
    peaks_odd = argrelmin(np.abs(psi_odd), order=order)[0]

    all_peaks = sorted(set(list(peaks_even) + list(peaks_odd)))
    return E_range[all_peaks]


# ============================================================
# Inverse power iteration (eigenvalue as fixed point)
# ============================================================

def inverse_power_iteration(H, sigma, n_iter=40, seed=42):
    """
    Inverse power iteration with shift sigma.

    (H - sigma*I)^{-1} applied iteratively converges to the
    eigenstate nearest sigma. This IS fixed-point iteration
    on the projective Hilbert space.

    Returns
    -------
    rayleigh_history : ndarray
        Rayleigh quotient at each iteration.
    """
    n = H.shape[0]
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(n)
    x /= np.linalg.norm(x)

    A = H - sigma * np.eye(n)
    rayleigh_history = np.zeros(n_iter)

    for i in range(n_iter):
        y = solve(A, x)
        x = y / np.linalg.norm(y)
        rayleigh_history[i] = x @ H @ x

    return rayleigh_history


def projective_power_iteration(H, x0, alpha=1.0, n_iter=200):
    """
    Alpha-transformed power iteration on projective Hilbert space.

        x_{n+1} = alpha * (H @ x_n / ||H @ x_n||) + (1-alpha) * x_n

    For alpha=1: standard power method (converges to largest |eigenvalue|).
    Other alpha values can select different eigenstates.

    Returns
    -------
    overlaps : ndarray, shape (n_iter, n_eigenstates)
        |<psi_k|x_n>|^2 at each iteration.
    eigenvalues : ndarray
        Eigenvalues of H.
    """
    eigenvalues, eigenvectors = eigh(H)
    n = len(eigenvalues)
    x = x0 / np.linalg.norm(x0)

    overlaps = np.zeros((n_iter, n))

    for i in range(n_iter):
        for k in range(n):
            overlaps[i, k] = np.abs(np.dot(eigenvectors[:, k], x))**2

        Hx = H @ x
        norm = np.linalg.norm(Hx)
        if norm < 1e-15:
            break
        Hx_normed = Hx / norm
        x_new = alpha * Hx_normed + (1 - alpha) * x
        x = x_new / np.linalg.norm(x_new)

    return overlaps, eigenvalues


def alpha_eigenstate_scan(H, x0, alphas, n_iter=300):
    """
    For each alpha, run projective power iteration and record
    which eigenstate the iteration converges to.

    Returns
    -------
    dominant : ndarray of int
        Index of dominant eigenstate for each alpha.
    quality : ndarray
        Overlap quality |<psi_k|x>|^2.
    eigenvalues : ndarray
    """
    eigenvalues, eigenvectors = eigh(H)
    n = len(eigenvalues)

    dominant = np.zeros(len(alphas), dtype=int)
    quality = np.zeros(len(alphas))

    for j, alpha in enumerate(alphas):
        x = x0 / np.linalg.norm(x0)
        for _ in range(n_iter):
            Hx = H @ x
            norm = np.linalg.norm(Hx)
            if norm < 1e-15:
                break
            Hx_normed = Hx / norm
            x_new = alpha * Hx_normed + (1 - alpha) * x
            x = x_new / np.linalg.norm(x_new)

        overlaps = np.array([np.abs(np.dot(eigenvectors[:, k], x))**2
                             for k in range(n)])
        dominant[j] = np.argmax(overlaps)
        quality[j] = overlaps[dominant[j]]

    return dominant, quality, eigenvalues
