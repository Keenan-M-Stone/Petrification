"""
Eigenstate-Fixed-Point Correspondence: Proof of Concept

Demonstrates three key ideas:
1. Eigenstates as fixed points of the projective power iteration
2. Turbiner's Riccati-Bloch nonlinearization as a discrete dynamical system
3. The alpha-transform as an eigenstate selector via stability control

The spectral bifurcation diagram: eigenvalues appear as bifurcation points
in the Riccati iteration parameter E.

Author: Keenan M Stone
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh


# ============================================================
# Part 1: Eigenstates as fixed points of projective iteration
# ============================================================

def projective_power_iteration(H, x0, alpha=1.0, n_iter=200):
    """
    Iterate x -> alpha * (H @ x / ||H @ x||) + (1 - alpha) * x, normalized.

    For alpha=1 this is the standard power method (converges to largest |eigenvalue|).
    For other alpha, it's the relaxed/alpha-transformed version that can
    converge to different eigenstates depending on alpha.

    Returns trajectory of overlap with each eigenstate.
    """
    eigenvalues, eigenvectors = eigh(H)
    n = len(eigenvalues)
    x = x0 / np.linalg.norm(x0)

    # Track overlap |<psi_k|x_n>|^2 with each eigenstate
    overlaps = np.zeros((n_iter, n))

    for i in range(n_iter):
        for k in range(n):
            overlaps[i, k] = np.abs(np.dot(eigenvectors[:, k], x))**2

        Hx = H @ x
        norm = np.linalg.norm(Hx)
        if norm < 1e-15:
            break
        Hx_normed = Hx / norm

        # Alpha-transform in projective space
        x_new = alpha * Hx_normed + (1 - alpha) * x
        x = x_new / np.linalg.norm(x_new)

    return overlaps, eigenvalues


# ============================================================
# Part 2: Riccati-Bloch discretization
# ============================================================

def riccati_iterate(V_func, E, x_grid, y0=0.0, alpha=1.0):
    """
    Integrate the Riccati-Bloch equation y' = y^2 - V(x) + E
    using forward Euler with alpha-relaxation.

    y = -psi'/psi (logarithmic derivative of the wavefunction).

    At eigenvalues, y(x) remains bounded.
    Away from eigenvalues, y(x) diverges.

    Returns y trajectory.
    """
    h = x_grid[1] - x_grid[0]
    y = np.copy(y0) if isinstance(y0, np.ndarray) else float(y0)
    trajectory = np.zeros(len(x_grid))
    trajectory[0] = y

    for i in range(1, len(x_grid)):
        x = x_grid[i - 1]
        # Riccati: y' = y^2 - 2V(x) + 2E   (from -psi''/2 + V*psi = E*psi)
        y_rhs = y**2 - 2.0 * V_func(x) + 2.0 * E
        y_euler = y + h * y_rhs

        # Alpha-transform: blend Euler step with identity
        y_new = alpha * y_euler + (1 - alpha) * y

        # Divergence clamp for visualization
        if abs(y_new) > 1e6:
            trajectory[i:] = np.sign(y_new) * 1e6
            break

        y = y_new
        trajectory[i] = y

    return trajectory


def spectral_bifurcation_diagram(V_func, E_range, x_grid, y0=0.0, alpha=1.0,
                                 tail_fraction=0.3):
    """
    For each trial energy E, run Riccati iteration and record
    the value of y at the tail end of the x-domain.

    At eigenvalues, y stays bounded -> visible as "bands".
    Away from eigenvalues, y diverges -> blank regions.

    This is the spectral bifurcation diagram.
    """
    n_tail = max(1, int(len(x_grid) * tail_fraction))
    results_E = []
    results_y = []

    for E in E_range:
        y_traj = riccati_iterate(V_func, E, x_grid, y0, alpha)
        tail = y_traj[-n_tail:]

        # Only keep bounded trajectories
        if np.all(np.abs(tail) < 1e5):
            for yval in tail[::max(1, len(tail) // 20)]:
                results_E.append(E)
                results_y.append(yval)

    return np.array(results_E), np.array(results_y)


# ============================================================
# Part 3: Alpha-scan — eigenstate selector
# ============================================================

def alpha_eigenstate_scan(H, x0, alphas, n_iter=300):
    """
    For each alpha, run projective power iteration and record
    which eigenstate the iteration converges to.

    Returns: for each alpha, the index of the dominant eigenstate.
    """
    eigenvalues, eigenvectors = eigh(H)
    n = len(eigenvalues)

    dominant = np.zeros(len(alphas), dtype=int)
    convergence_quality = np.zeros(len(alphas))

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

        # Which eigenstate did we land on?
        overlaps = np.array([np.abs(np.dot(eigenvectors[:, k], x))**2
                             for k in range(n)])
        dominant[j] = np.argmax(overlaps)
        convergence_quality[j] = overlaps[dominant[j]]

    return dominant, convergence_quality, eigenvalues


# ============================================================
# Part 4: Potentials
# ============================================================

def harmonic_potential(x, omega=1.0):
    return 0.5 * omega**2 * x**2


def anharmonic_potential(x, g=0.1):
    """Quartic anharmonic oscillator: V = x^2 + g*x^4"""
    return x**2 + g * x**4


def double_well(x, g=1.0):
    """Double well: V = x^2(1 - g*x)^2"""
    return x**2 * (1 - g * x)**2


def finite_well(x, V0=50.0, a=1.0):
    """Finite square well."""
    return np.where(np.abs(x) < a, 0.0, V0)


def discretize_hamiltonian(V_func, x_grid):
    """
    Build the Hamiltonian matrix using finite differences.
    H = -d^2/dx^2 + V(x) in natural units (hbar=m=1).
    """
    n = len(x_grid)
    h = x_grid[1] - x_grid[0]

    # Kinetic energy: -0.5 * d^2/dx^2 via finite differences (hbar = m = 1)
    T = np.zeros((n, n))
    for i in range(n):
        T[i, i] = 1.0 / h**2
        if i > 0:
            T[i, i - 1] = -0.5 / h**2
        if i < n - 1:
            T[i, i + 1] = -0.5 / h**2

    # Potential energy
    V = np.diag(V_func(x_grid))

    return T + V


# ============================================================
# Driver: generate all four key plots
# ============================================================

def main():
    print("=" * 70)
    print("EIGENSTATE-FIXED-POINT CORRESPONDENCE")
    print("Spectral Bifurcation via Riccati-Bloch Dynamics")
    print("=" * 70)

    # --- Setup ---
    N = 500
    x_grid = np.linspace(-8, 8, N)
    V_func = harmonic_potential

    H = discretize_hamiltonian(V_func, x_grid)
    eigenvalues, eigenvectors = eigh(H)

    print(f"\nFirst 8 eigenvalues (harmonic oscillator, N={N}):")
    exact = np.arange(8) + 0.5
    for i in range(min(8, len(eigenvalues))):
        print(f"  E_{i} = {eigenvalues[i]:.6f}  (exact: {exact[i]:.1f})")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle("Eigenstate–Fixed Point Correspondence", fontsize=15, y=0.98)

    # ---- Plot 1: Projective power iteration convergence ----
    ax = axes[0, 0]
    np.random.seed(42)
    x0 = np.random.randn(N)
    overlaps, evals = projective_power_iteration(H, x0, alpha=1.0, n_iter=60)
    for k in range(min(5, overlaps.shape[1])):
        ax.plot(overlaps[:, k], label=f"$|\\langle\\psi_{k}|x_n\\rangle|^2$")
    ax.set_xlabel("Iteration $n$")
    ax.set_ylabel("Overlap with eigenstate")
    ax.set_title("Power iteration converges to\nfixed point = dominant eigenstate")
    ax.legend(fontsize=8)
    ax.set_ylim(-0.05, 1.05)

    # ---- Plot 2: Spectral bifurcation diagram ----
    ax = axes[0, 1]
    E_range = np.linspace(0.0, 8.5, 3000)
    x_riccati = np.linspace(0.01, 5.5, 600)

    bif_E, bif_y = spectral_bifurcation_diagram(
        V_func, E_range, x_riccati, y0=0.01, alpha=1.0, tail_fraction=0.2
    )

    if len(bif_E) > 0:
        ax.scatter(bif_E, bif_y, s=0.1, c='steelblue', alpha=0.3)

    # Mark exact eigenvalues
    for i, ev in enumerate(eigenvalues[:8]):
        if 0 <= ev <= 8.5:
            ax.axvline(ev, color='red', alpha=0.4, linewidth=0.8,
                       label=f"$E_{i}={ev:.2f}$" if i < 4 else None)

    ax.set_xlabel("Trial energy $E$")
    ax.set_ylabel("$y(x)$ at tail")
    ax.set_title("Spectral bifurcation diagram\n(Riccati-Bloch, harmonic oscillator)")
    ax.set_ylim(-10, 10)
    ax.legend(fontsize=7, loc="upper left")

    # ---- Plot 3: Alpha-transform eigenstate selection ----
    ax = axes[1, 0]
    # Use a smaller system for clarity
    N_small = 60
    x_small = np.linspace(-5, 5, N_small)
    H_small = discretize_hamiltonian(V_func, x_small)
    evals_small, evecs_small = eigh(H_small)

    alphas = np.linspace(-2.0, 2.0, 500)
    np.random.seed(7)
    x0_small = np.random.randn(N_small)
    dom, quality, ev_small = alpha_eigenstate_scan(H_small, x0_small, alphas, n_iter=500)

    scatter = ax.scatter(alphas, dom, c=quality, cmap='viridis', s=3, vmin=0, vmax=1)
    ax.set_xlabel("$\\alpha$")
    ax.set_ylabel("Dominant eigenstate index")
    ax.set_title("$\\alpha$-transform selects different\neigenstates (color = convergence)")
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("$|\\langle\\psi_k|x\\rangle|^2$")

    # ---- Plot 4: Riccati trajectories at/near eigenvalue ----
    ax = axes[1, 1]
    E_exact = eigenvalues[0]  # ground state

    offsets = [-0.3, -0.1, 0.0, 0.1, 0.3]
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#ff7f0e', '#d62728']
    styles = ['--', '-.', '-', '-.', '--']

    for offset, color, style in zip(offsets, colors, styles):
        E_trial = E_exact + offset
        y_traj = riccati_iterate(V_func, E_trial, x_riccati, y0=0.01, alpha=1.0)
        label = f"$E = E_0 {'+' if offset >= 0 else ''}{offset:.1f}$"
        ax.plot(x_riccati, y_traj, color=color, linestyle=style,
                linewidth=1.5 if offset == 0 else 1.0, label=label)

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y(x) = -\\psi'/\\psi$")
    ax.set_title(f"Riccati trajectories near $E_0={E_exact:.3f}$\n"
                 "Bounded iff $E$ = eigenvalue")
    ax.set_ylim(-15, 15)
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("/home/lemma137/dev/Petrification/eigen_fixedpoint_demo.png",
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nPlot saved to Petrification/eigen_fixedpoint_demo.png")

    # ---- Numerical verification ----
    print("\n" + "=" * 70)
    print("VERIFICATION: Riccati divergence rate encodes spectral gap")
    print("=" * 70)

    E_test = np.linspace(0.0, 5.0, 500)
    max_y = []
    for E in E_test:
        y_traj = riccati_iterate(V_func, E, x_riccati, y0=0.01, alpha=1.0)
        max_y.append(np.max(np.abs(y_traj)))

    max_y = np.array(max_y)

    # Find minima of max|y| — these should be near eigenvalues
    from scipy.signal import argrelmin
    local_mins = argrelmin(max_y, order=10)[0]
    print("\nRiccati-detected eigenvalues (minima of max|y|):")
    for idx in local_mins:
        E_detected = E_test[idx]
        # Find closest exact eigenvalue
        closest = eigenvalues[np.argmin(np.abs(eigenvalues - E_detected))]
        error = abs(E_detected - closest)
        print(f"  E_detected = {E_detected:.4f}, closest exact = {closest:.4f}, "
              f"error = {error:.4f}")


if __name__ == "__main__":
    main()
