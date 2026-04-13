"""
Eigenstate-Fixed-Point Correspondence: Proof of Concept (v2)

Demonstrates three key ideas:
1. Eigenstates as fixed points: inverse power iteration with shift σ converges
   to the eigenstate nearest σ — a fixed point of (H - σI)^{-1}.
2. Spectral bifurcation diagram: the Riccati-Bloch equation y' = y² - 2V + 2E
   produces bounded solutions ONLY at eigenvalues E_n. The "boundedness"
   as a function of E creates a bifurcation-like diagram.
3. Alpha-transform as eigenstate selector: parametrized relaxation in the
   Riccati iteration can stabilize different eigenstates.

Author: Keenan M Stone
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh, solve
from scipy.integrate import solve_ivp


# ============================================================
# Potentials
# ============================================================

def harmonic_potential(x, omega=1.0):
    return 0.5 * omega**2 * x**2


def anharmonic_potential(x, g=0.1):
    return x**2 + g * x**4


def double_well(x, a=1.0):
    return (x**2 - a**2)**2


# ============================================================
# Discretized Hamiltonian
# ============================================================

def discretize_hamiltonian(V_func, x_grid):
    n = len(x_grid)
    h = x_grid[1] - x_grid[0]
    T = np.zeros((n, n))
    for i in range(n):
        T[i, i] = 1.0 / h**2
        if i > 0:
            T[i, i - 1] = -0.5 / h**2
        if i < n - 1:
            T[i, i + 1] = -0.5 / h**2
    V = np.diag(V_func(x_grid))
    return T + V


# ============================================================
# Demo 1: Inverse power method = fixed-point iteration
# ============================================================

def inverse_power_demo(H, eigenvalues, eigenvectors, sigma, n_iter=40):
    """
    Run inverse power iteration with shift sigma.
    (H - sigma*I)^{-1} x_n / ||(H - sigma*I)^{-1} x_n|| -> eigenstate nearest sigma.
    This IS fixed-point iteration on the projective space.
    """
    n = H.shape[0]
    np.random.seed(42)
    x = np.random.randn(n)
    x /= np.linalg.norm(x)

    A = H - sigma * np.eye(n)

    rayleigh_history = []
    for _ in range(n_iter):
        y = solve(A, x)
        x = y / np.linalg.norm(y)
        rq = x @ H @ x
        rayleigh_history.append(rq)

    return np.array(rayleigh_history)


# ============================================================
# Demo 2: Riccati-Bloch via high-order integration
# ============================================================

def riccati_rhs(x, y, V_func, E):
    """y' = y^2 - 2V(x) + 2E"""
    return y**2 - 2.0 * V_func(x) + 2.0 * E


def riccati_solve(V_func, E, x_span, y0=0.0, max_y=100.0):
    """
    Integrate the Riccati-Bloch equation y' = y^2 - 2V + 2E from
    x_span[0] to x_span[1] using RK45. Return trajectory.

    Stops if |y| > max_y (divergence).
    """
    def rhs(x, y):
        return [riccati_rhs(x, y[0], V_func, E)]

    def divergence_event(x, y):
        return max_y - abs(y[0])
    divergence_event.terminal = True
    divergence_event.direction = -1

    sol = solve_ivp(rhs, x_span, [y0], method='RK45',
                    events=[divergence_event],
                    max_step=0.01, dense_output=True,
                    rtol=1e-8, atol=1e-10)
    return sol


def spectral_scan(V_func, E_range, x_span=(0.01, 8.0)):
    """
    For each trial E, integrate the Schrödinger equation forward as a
    shooting problem and record psi(x_max).

    At eigenvalues: psi(x_max) -> 0 (boundary condition satisfied).
    Away from eigenvalues: psi(x_max) diverges.

    This is a fixed-point problem: find E such that the shooting map
    S(E) = psi(x_max; E) has a fixed point S(E) = 0.

    We plot log10(|psi(x_max)|) to visualize the sharp dips.

    For the Riccati connection: psi = exp(-int y dx), so
    psi(x_max) -> 0  iff  int y dx -> +infinity in a controlled way.
    """
    results = np.zeros(len(E_range))

    for j, E in enumerate(E_range):
        def schrodinger(x, state):
            psi, dpsi = state
            # -psi''/2 + V*psi = E*psi  =>  psi'' = 2*(V-E)*psi
            return [dpsi, 2.0 * (V_func(x) - E) * psi]

        # Even-parity initial condition: psi(0)=1, psi'(0)=0
        sol = solve_ivp(schrodinger, [x_span[0], x_span[1]], [1.0, 0.0],
                        method='RK45', max_step=0.02, rtol=1e-10, atol=1e-12,
                        dense_output=True)
        psi_end = sol.y[0, -1]
        results[j] = psi_end

    return results


def spectral_scan_odd(V_func, E_range, x_span=(0.01, 8.0)):
    """Same as spectral_scan but for odd-parity states: psi(0)=0, psi'(0)=1."""
    results = np.zeros(len(E_range))

    for j, E in enumerate(E_range):
        def schrodinger(x, state):
            psi, dpsi = state
            return [dpsi, 2.0 * (V_func(x) - E) * psi]

        sol = solve_ivp(schrodinger, [x_span[0], x_span[1]], [0.0, 1.0],
                        method='RK45', max_step=0.02, rtol=1e-10, atol=1e-12,
                        dense_output=True)
        psi_end = sol.y[0, -1]
        results[j] = psi_end

    return results


# ============================================================
# Demo 3: Alpha-relaxed Riccati
# ============================================================

def riccati_alpha(V_func, E, x_grid, y0=0.0, alpha=1.0):
    """
    Integrate Riccati with alpha-relaxation using Euler for explicit control.
    y_{n+1} = alpha * [y_n + h*(y_n^2 - 2V + 2E)] + (1-alpha) * y_n
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

        if abs(y_new) > 1e4:
            trajectory[i:] = np.sign(y_new) * 1e4
            break
        y = y_new
        trajectory[i] = y

    return trajectory


# ============================================================
# Driver
# ============================================================

def main():
    print("=" * 70)
    print("EIGENSTATE-FIXED-POINT CORRESPONDENCE (v2)")
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
    for i in range(8):
        print(f"  E_{i} = {eigenvalues[i]:.6f}  (exact: {exact[i]:.1f})")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle("Eigenstate–Fixed Point Correspondence", fontsize=15, y=0.98)

    # ====================================================
    # Plot 1: Inverse power method as fixed-point iteration
    # ====================================================
    ax = axes[0, 0]
    shifts = [0.4, 1.6, 2.4, 3.6, 4.6]
    colors_shifts = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    for sigma, color in zip(shifts, colors_shifts):
        rq_hist = inverse_power_demo(H, eigenvalues, eigenvectors, sigma)
        nearest_idx = np.argmin(np.abs(eigenvalues - sigma))
        nearest_E = eigenvalues[nearest_idx]
        ax.plot(rq_hist, color=color, linewidth=1.8,
                label=f"$\\sigma={sigma:.1f} \\to E_{nearest_idx}={nearest_E:.2f}$")

    for i in range(5):
        ax.axhline(eigenvalues[i], color='gray', alpha=0.3, linewidth=0.8)

    ax.set_xlabel("Iteration $n$", fontsize=11)
    ax.set_ylabel("Rayleigh quotient $\\langle x|H|x\\rangle$", fontsize=11)
    ax.set_title("Inverse power iteration = fixed-point\nconvergence to nearest eigenstate")
    ax.legend(fontsize=8, loc='right')
    ax.set_ylim(0, 5.5)

    # ====================================================
    # Plot 2: Spectral scan — shooting method
    # ====================================================
    ax = axes[0, 1]
    E_range = np.linspace(0.01, 8.5, 2000)
    x_span = (0.0, 7.0)

    psi_even = spectral_scan(V_func, E_range, x_span)
    psi_odd = spectral_scan_odd(V_func, E_range, x_span)

    # Plot log10(|psi|) — dips are eigenvalues
    log_even = np.log10(np.abs(psi_even) + 1e-20)
    log_odd = np.log10(np.abs(psi_odd) + 1e-20)

    ax.plot(E_range, log_even, 'steelblue', linewidth=1.0, label='Even parity')
    ax.plot(E_range, log_odd, '#ff7f0e', linewidth=1.0, label='Odd parity')

    # Mark eigenvalues
    for i in range(8):
        ev = eigenvalues[i]
        if ev <= 8.5:
            ax.axvline(ev, color='red', alpha=0.4, linewidth=0.8,
                       linestyle='--')

    ax.set_xlabel("Trial energy $E$", fontsize=11)
    ax.set_ylabel("$\\log_{10}|\\psi(x_{\\max}, E)|$", fontsize=11)
    ax.set_title("Spectral shooting: $\\psi(x_{\\max})\\to 0$\nonly at eigenvalues (= fixed points of $S(E)$)")
    ax.legend(fontsize=9)
    ax.set_ylim(-18, 20)

    # ====================================================
    # Plot 3: Riccati trajectories at/near eigenvalue
    # ====================================================
    ax = axes[1, 0]
    E_exact = eigenvalues[0]
    x_riccati_grid = np.linspace(0.01, 6.0, 2000)

    offsets = [-0.3, -0.1, -0.01, 0.0, 0.01, 0.1, 0.3]
    colors_traj = ['#d62728', '#ff7f0e', '#bcbd22', '#2ca02c', '#bcbd22', '#ff7f0e', '#d62728']
    widths = [0.8, 1.0, 1.2, 2.5, 1.2, 1.0, 0.8]

    for offset, color, w in zip(offsets, colors_traj, widths):
        E_trial = E_exact + offset
        label = f"$E = E_0 {'+' if offset >= 0 else ''}{offset}$"
        sol = riccati_solve(V_func, E_trial, (0.01, 6.0), y0=0.01)
        ax.plot(sol.t, sol.y[0], color=color, linewidth=w, label=label)

    ax.set_xlabel("$x$", fontsize=11)
    ax.set_ylabel("$y(x) = -\\psi'/\\psi$", fontsize=11)
    ax.set_title(f"Riccati trajectories near $E_0={E_exact:.4f}$\n"
                 "Bounded (green) = eigenvalue (fixed point)")
    ax.set_ylim(-20, 20)
    ax.legend(fontsize=7, loc='lower left')

    # ====================================================
    # Plot 4: Alpha-relaxation of Riccati at E=E_1
    # ====================================================
    ax = axes[1, 1]
    E1 = eigenvalues[1]  # first excited state
    x_alpha_grid = np.linspace(0.01, 5.0, 1500)

    alphas_demo = [0.3, 0.5, 0.7, 1.0, 1.3]
    colors_alpha = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    for alpha_val, color in zip(alphas_demo, colors_alpha):
        y_traj = riccati_alpha(V_func, E1, x_alpha_grid, y0=0.01, alpha=alpha_val)
        ax.plot(x_alpha_grid, y_traj, color=color, linewidth=1.3,
                label=f"$\\alpha = {alpha_val}$")

    ax.set_xlabel("$x$", fontsize=11)
    ax.set_ylabel("$y(x)$", fontsize=11)
    ax.set_title(f"Alpha-relaxed Riccati at $E_1 = {E1:.4f}$\n"
                 "$\\alpha$ controls stability of the trajectory")
    ax.set_ylim(-20, 20)
    ax.legend(fontsize=9)

    plt.tight_layout()
    outpath = "/home/lemma137/dev/Petrification/eigen_fixedpoint_demo.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved to {outpath}")

    # ====================================================
    # Quantitative check: eigenvalues from Riccati peaks
    # ====================================================
    print("\n" + "=" * 70)
    print("EIGENVALUE DETECTION VIA RICCATI PENETRATION DEPTH")
    print("=" * 70)

    E_fine = np.linspace(0.01, 8.0, 5000)
    psi_e = spectral_scan(V_func, E_fine, (0.0, 7.0))
    psi_o = spectral_scan_odd(V_func, E_fine, (0.0, 7.0))

    # Eigenvalues = sign changes (zeros of psi_end)
    from scipy.signal import argrelmin
    # Even-parity eigenvalues: dips in |psi_even|
    log_e = np.log10(np.abs(psi_e) + 1e-20)
    log_o = np.log10(np.abs(psi_o) + 1e-20)
    peaks_e = argrelmin(np.abs(psi_e), order=20)[0]
    peaks_o = argrelmin(np.abs(psi_o), order=20)[0]
    peaks = sorted(list(peaks_e) + list(peaks_o), key=lambda i: E_fine[i])
    print(f"\nDetected {len(peaks)} eigenvalues from Riccati penetration peaks:")
    for idx in peaks:
        E_det = E_fine[idx]
        near = np.argmin(np.abs(eigenvalues - E_det))
        err = abs(E_det - eigenvalues[near])
        print(f"  E_detected = {E_det:.5f}, "
              f"E_exact = {eigenvalues[near]:.5f}, "
              f"|error| = {err:.5f}")

    # ====================================================
    # Anharmonic oscillator test
    # ====================================================
    print("\n" + "=" * 70)
    print("ANHARMONIC OSCILLATOR: V = x^2 + 0.1*x^4")
    print("=" * 70)

    V_anh = lambda x: anharmonic_potential(x, g=0.1)
    H_anh = discretize_hamiltonian(V_anh, x_grid)
    evals_anh, _ = eigh(H_anh)

    print("\nFirst 6 eigenvalues (from matrix diagonalization):")
    for i in range(6):
        print(f"  E_{i} = {evals_anh[i]:.6f}")

    E_anh_range = np.linspace(0.01, 12.0, 3000)
    psi_anh_e = spectral_scan(V_anh, E_anh_range, (0.0, 7.0))
    psi_anh_o = spectral_scan_odd(V_anh, E_anh_range, (0.0, 7.0))

    peaks_ae = argrelmin(np.abs(psi_anh_e), order=10)[0]
    peaks_ao = argrelmin(np.abs(psi_anh_o), order=10)[0]
    peaks_anh = sorted(list(peaks_ae) + list(peaks_ao),
                       key=lambda i: E_anh_range[i])
    print(f"\nDetected {len(peaks_anh)} eigenvalues from Riccati scan:")
    for idx in peaks_anh:
        E_det = E_anh_range[idx]
        near = np.argmin(np.abs(evals_anh - E_det))
        err = abs(E_det - evals_anh[near])
        print(f"  E_detected = {E_det:.4f}, "
              f"E_exact = {evals_anh[near]:.4f}, "
              f"|error| = {err:.4f}")

    # ---- Anharmonic bifurcation plot ----
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    log_ae = np.log10(np.abs(psi_anh_e) + 1e-20)
    log_ao = np.log10(np.abs(psi_anh_o) + 1e-20)
    ax2.plot(E_anh_range, log_ae, 'steelblue', linewidth=0.8, label='Even')
    ax2.plot(E_anh_range, log_ao, '#ff7f0e', linewidth=0.8, label='Odd')
    for i in range(6):
        ev = evals_anh[i]
        if ev <= 12.0:
            ax2.axvline(ev, color='red', alpha=0.5, linewidth=0.8, linestyle='--')
    ax2.set_xlabel("Trial energy $E$", fontsize=12)
    ax2.set_ylabel("$\\log_{10}|\\psi(x_{\\max})|$", fontsize=12)
    ax2.set_title("Spectral shooting: anharmonic oscillator $V = x^2 + 0.1x^4$",
                  fontsize=13)
    ax2.legend()
    ax2.set_ylim(-18, 20)
    outpath2 = "/home/lemma137/dev/Petrification/anharmonic_bifurcation.png"
    plt.tight_layout()
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nAnharmonic plot saved to {outpath2}")


if __name__ == "__main__":
    main()
