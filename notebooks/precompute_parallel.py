"""
Parallel precomputation of expensive notebook results using Dask.

Computes and caches results for:
  1. lyapunov_alpha_relationship.ipynb  (Lyapunov phase diagram + period map)
  2. piecewise_alpha_dyson.ipynb        (stability landscapes + piecewise sweeps)
  3. rp_resonance_basis.ipynb           (Ulam + Chebyshev + Jacobi transfer operators)

Results saved as .npz files in notebooks/cache/
"""

import sys
sys.path.insert(0, '..')

import os
import time
import warnings
import numpy as np
from scipy.special import eval_chebyt, eval_jacobi
from scipy import integrate

from petrification.maps import logistic
from petrification.iteration import iterate, iterate_transformed
from petrification.quantum import frobenius_perron_matrix
from petrification.dyson import (
    scalar_dyson_iterate, scalar_dyson_exact, scalar_dyson_stability,
    spectral_dyson_scan, quadratic_self_energy, cubic_self_energy
)

CACHE_DIR = os.path.join(os.path.dirname(__file__), 'cache')
os.makedirs(CACHE_DIR, exist_ok=True)


# ================================================================
# Worker functions (must be picklable — no closures)
# ================================================================

def _lyapunov_single(args):
    """Compute Lyapunov exponent at a single (alpha, a) point."""
    alpha_val, a_val, x0, n_iter, n_discard = args
    try:
        x = float(x0)
        traj = np.empty(n_iter + 1)
        traj[0] = x
        for k in range(n_iter):
            fx = a_val * x * (1 - x)
            if alpha_val is not None:
                x = alpha_val * fx + (1 - alpha_val) * x
            else:
                x = fx
            # Early bailout on divergence
            if abs(x) > 1e10 or x != x:  # x != x is fast NaN check
                return np.nan
            traj[k + 1] = x

        traj_use = traj[n_discard:]
        if len(traj_use) < 100:
            return np.nan

        fp = a_val * (1 - 2 * traj_use[:-1])
        if alpha_val is not None and alpha_val != 1.0:
            gp = 1 + alpha_val * (fp - 1)
        else:
            gp = fp
        gp_abs = np.maximum(np.abs(gp), 1e-30)
        return float(np.mean(np.log(gp_abs)))
    except (OverflowError, FloatingPointError, ValueError):
        return np.nan


def _period_single(args):
    """Detect period at a single (alpha, a) point."""
    alpha_val, a_val, x0, n_iter, tol, max_period = args
    try:
        x = float(x0)
        traj = np.empty(n_iter + 1)
        traj[0] = x
        for k in range(n_iter):
            fx = a_val * x * (1 - x)
            x = alpha_val * fx + (1 - alpha_val) * x
            if abs(x) > 1e10 or x != x:
                return -1
            traj[k + 1] = x

        tail = traj[-500:]
        x_final = tail[-1]
        for p in range(1, max_period + 1):
            if abs(tail[-(p + 1)] - x_final) < tol:
                if all(abs(tail[-(p + 1 + k)] - tail[-(1 + k)]) < tol
                       for k in range(min(p, 10))):
                    return p
        return 0  # chaotic
    except (OverflowError, FloatingPointError, ValueError):
        return -1


def _ulam_eigenvalues(N_val):
    """Compute Ulam transfer operator eigenvalues at resolution N."""
    L, bins = frobenius_perron_matrix(logistic, 4.0, N=N_val)
    eigs = np.linalg.eigvals(L)
    return np.sort(np.abs(eigs))[::-1][:6]


def _chebyshev_transfer_weighted(N, n_quad, a_param=4.0):
    """Compute transfer operator in Chebyshev basis."""
    k = np.arange(1, n_quad + 1)
    t_nodes = np.cos((2 * k - 1) * np.pi / (2 * n_quad))
    x_nodes = (t_nodes + 1) / 2
    w_nodes = np.ones(n_quad) / n_quad

    basis = np.zeros((N, n_quad))
    for n in range(N):
        basis[n, :] = eval_chebyt(n, t_nodes)

    L_basis = np.zeros((N, n_quad))
    a = a_param
    for j, x in enumerate(x_nodes):
        disc = 1.0 - x
        if disc < 0:
            continue
        sqrt_disc = np.sqrt(disc)
        preimages = [(1 + sqrt_disc) / 2, (1 - sqrt_disc) / 2]
        for y in preimages:
            if y < 0 or y > 1:
                continue
            fp = abs(a * (1 - 2 * y))
            if fp < 1e-15:
                continue
            t_y = 2 * y - 1
            for n in range(N):
                L_basis[n, j] += eval_chebyt(n, t_y) / fp

    L_mat = np.zeros((N, N))
    for m in range(N):
        for n in range(N):
            L_mat[m, n] = np.sum(w_nodes * basis[m, :] * L_basis[n, :])
    return L_mat


def _jacobi_transfer_weighted(args):
    """Compute transfer operator in Jacobi basis."""
    N, alpha_j, beta_j, n_quad, a_param = args
    from numpy.polynomial.legendre import leggauss
    gl_nodes, gl_weights = leggauss(n_quad)
    x_nodes = (gl_nodes + 1) / 2

    jacobi_weight = np.zeros(n_quad)
    for j in range(n_quad):
        x = x_nodes[j]
        if x > 1e-12 and x < 1 - 1e-12:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                w = x**alpha_j * (1 - x)**beta_j
                if np.isfinite(w):
                    jacobi_weight[j] = w

    w_eff = gl_weights * jacobi_weight / 2
    norm_total = np.sum(w_eff)
    if norm_total > 0:
        w_eff = w_eff / norm_total

    basis = np.zeros((N, n_quad))
    for n in range(N):
        basis[n, :] = eval_jacobi(n, alpha_j, beta_j, 2 * x_nodes - 1)

    L_basis = np.zeros((N, n_quad))
    a = a_param
    for j, x in enumerate(x_nodes):
        disc = 1.0 - x
        if disc < 0:
            continue
        sqrt_disc = np.sqrt(disc)
        preimages = [(1 + sqrt_disc) / 2, (1 - sqrt_disc) / 2]
        for y in preimages:
            if y < 1e-15 or y > 1 - 1e-15:
                continue
            fp = abs(a * (1 - 2 * y))
            if fp < 1e-15:
                continue
            t_y = 2 * y - 1
            for n in range(N):
                L_basis[n, j] += eval_jacobi(n, alpha_j, beta_j, t_y) / fp

    L_mat = np.zeros((N, N))
    for m in range(N):
        for n in range(N):
            L_mat[m, n] = np.sum(w_eff * basis[m, :] * L_basis[n, :])
    return L_mat


def _make_sigma(sigma_type, U, V):
    """Create a self-energy function (top-level for pickling)."""
    if sigma_type == 'cubic':
        return cubic_self_energy(U, V)
    elif sigma_type == 'multipole':
        return _multipole_self_energy(U, V)
    else:
        return quadratic_self_energy(U)


def _multipole_self_energy(U, V):
    """Σ(G) = U²G/(1 - V²G²) — top-level function for pickling."""
    def sigma(G):
        denom = 1 - V**2 * G**2
        if abs(denom) < 1e-20:
            return U**2 * G
        return U**2 * G / denom
    return sigma


def _dyson_stability_single(args):
    """Numerical stability for arbitrary self-energy."""
    omega_real, eps0, sigma_type, U, V = args
    omega = complex(omega_real, 0.05)

    if sigma_type == 'quadratic':
        s, _, a_opt = scalar_dyson_stability(omega, eps0, U)
        return float(s), float(np.real(a_opt))

    sigma_func = _make_sigma(sigma_type, U, V)
    try:
        G0 = 1.0 / (omega - eps0)
        hist, conv = scalar_dyson_iterate(omega, eps0, sigma_func, G0,
                                          alpha=0.5, n_iter=500, tol=1e-14)
        G_star = hist[-1]

        def F(G):
            return 1.0 / (omega - eps0 - sigma_func(G))

        dG = 1e-8 * (1 + 1j)
        F_prime = (F(G_star + dG) - F(G_star - dG)) / (2 * dG)
        stab = abs(F_prime)
        a_opt = np.real(1.0 / (1.0 - F_prime)) if abs(1.0 - F_prime) > 1e-10 else 1.0
        return float(stab), float(np.clip(a_opt, 0, 1))
    except Exception:
        return float('nan'), 0.5


def _piecewise_window_sweep(args):
    """Find optimal alpha for a single frequency window."""
    omega_slice_real, eps0, sigma_type, U, V, alpha_candidates = args
    eta = 0.05
    sigma_func = _make_sigma(sigma_type, U, V)

    best_a, best_score = 0.5, 1e10
    for a_try in alpha_candidates:
        total_iters, n_conv = 0, 0
        for w_r in omega_slice_real:
            omega = complex(w_r, eta)
            G0 = 1.0 / (omega - eps0)
            hist, conv = scalar_dyson_iterate(omega, eps0, sigma_func, G0,
                                              alpha=a_try, n_iter=500, tol=1e-12)
            total_iters += len(hist) - 1
            n_conv += int(conv)
        score = total_iters / max(n_conv, 1) + 500 * (len(omega_slice_real) - n_conv)
        if score < best_score:
            best_a, best_score = a_try, score
    return best_a


# ================================================================
# Main precomputation
# ================================================================

def precompute_all():
    from dask.distributed import Client, LocalCluster

    # Suppress expected warnings from chaotic trajectory edge cases
    warnings.filterwarnings('ignore', category=RuntimeWarning,
                            message='.*overflow.*|.*invalid value.*|.*divide by zero.*')

    cluster = LocalCluster(n_workers=4, threads_per_worker=1, memory_limit='2GB')
    client = Client(cluster)
    print(f"Dask dashboard: {client.dashboard_link}")
    print(f"Workers: {len(client.scheduler_info()['workers'])}")

    t0 = time.time()

    # ============================================================
    # 1. LYAPUNOV NOTEBOOK
    # ============================================================
    print("\n" + "=" * 60)
    print("1. Lyapunov notebook computations")
    print("=" * 60)

    # §3: Regular vs inverted scan (500 points each)
    print("  §3: Regular vs inverted scan...")
    a_values = np.linspace(2.5, 4.0, 500)

    args_reg = [(None, a, 0.4, 10000, 1000) for a in a_values]
    args_inv = [(-1.0, a, 0.4, 10000, 1000) for a in a_values]

    futures_reg = client.map(_lyapunov_single, args_reg)
    futures_inv = client.map(_lyapunov_single, args_inv)
    lyap_regular = np.array(client.gather(futures_reg))
    lyap_inverted = np.array(client.gather(futures_inv))

    # §4: Continuous alpha sweep at 4 parameter values
    print("  §4: Alpha sweep at 4 parameter values...")
    alpha_range = np.linspace(-2.0, 2.0, 200)
    a_test = [3.2, 3.57, 3.83, 4.0]
    lyap_alpha_sweep = {}

    for a in a_test:
        args_sweep = [(alpha, a, 0.4, 12000, 2000) for alpha in alpha_range]
        futures = client.map(_lyapunov_single, args_sweep)
        lyap_alpha_sweep[a] = np.array(client.gather(futures))

    # §5: Phase diagram (200 × 300 = 60,000 points) — THE BIG ONE
    print("  §5: Phase diagram (200×300 = 60,000 points)...")
    n_alpha = 200
    n_a = 300
    alpha_range_2d = np.linspace(-1.5, 1.5, n_alpha)
    a_range_2d = np.linspace(2.5, 4.0, n_a)

    args_2d = [(alpha, a, 0.4, 6000, 2000)
               for alpha in alpha_range_2d for a in a_range_2d]
    futures_2d = client.map(_lyapunov_single, args_2d, batch_size=500)
    results_2d = client.gather(futures_2d)
    lyap_map = np.array(results_2d).reshape(n_alpha, n_a)

    # §6: Theoretical prediction at a=4 + computed comparison
    print("  §6: Theory vs computation at a=4...")
    alpha_test_theory = np.linspace(-1.5, 2.0, 50)
    args_theory = [(a if a != 1.0 else None, 4.0, 0.4, 50000, 1000) for a in alpha_test_theory]
    futures_theory = client.map(_lyapunov_single, args_theory)
    lyap_computed_theory = np.array(client.gather(futures_theory))

    # §7: Period map (150 × 200 = 30,000 points)
    print("  §7: Period map (150×200 = 30,000 points)...")
    n_alpha_p = 150
    n_a_p = 200
    alpha_range_p = np.linspace(-1.5, 1.5, n_alpha_p)
    a_range_p = np.linspace(2.5, 4.0, n_a_p)

    args_period = [(alpha, a, 0.4, 5000, 1e-6, 256)
                   for alpha in alpha_range_p for a in a_range_p]
    futures_period = client.map(_period_single, args_period, batch_size=500)
    results_period = client.gather(futures_period)
    period_map = np.array(results_period).reshape(n_alpha_p, n_a_p)

    np.savez(os.path.join(CACHE_DIR, 'lyapunov_results.npz'),
             a_values=a_values, lyap_regular=lyap_regular, lyap_inverted=lyap_inverted,
             alpha_range=alpha_range, a_test=np.array(a_test),
             lyap_sweep_3_2=lyap_alpha_sweep[3.2], lyap_sweep_3_57=lyap_alpha_sweep[3.57],
             lyap_sweep_3_83=lyap_alpha_sweep[3.83], lyap_sweep_4_0=lyap_alpha_sweep[4.0],
             alpha_range_2d=alpha_range_2d, a_range_2d=a_range_2d, lyap_map=lyap_map,
             alpha_test_theory=alpha_test_theory, lyap_computed_theory=lyap_computed_theory,
             alpha_range_p=alpha_range_p, a_range_p=a_range_p, period_map=period_map)
    print("  -> Saved lyapunov_results.npz")

    # ============================================================
    # 2. PIECEWISE DYSON NOTEBOOK
    # ============================================================
    print("\n" + "=" * 60)
    print("2. Piecewise Dyson notebook computations")
    print("=" * 60)

    eta = 0.05
    eps0 = 0.0
    U_pw = 3.0
    omega_grid_500 = np.linspace(-5, 5, 500)
    omega_grid_200 = np.linspace(-5, 5, 200)

    # §2: Stability landscapes for 3 self-energies
    print("  §2: Stability landscapes...")
    for stype, V in [('quadratic', 0), ('cubic', 1.0), ('multipole', 0.8)]:
        args_stab = [(w, eps0, stype, U_pw, V) for w in omega_grid_500]
        futures_stab = client.map(_dyson_stability_single, args_stab)
        results_stab = client.gather(futures_stab)
        stab_arr = np.array([r[0] for r in results_stab])
        aopt_arr = np.array([r[1] for r in results_stab])
        np.savez(os.path.join(CACHE_DIR, f'stability_{stype}.npz'),
                 omega=omega_grid_500, stability=stab_arr, alpha_opt=aopt_arr)
        print(f"    -> Saved stability_{stype}.npz")

    # §3-5: Piecewise window optimization for cubic & multipole
    print("  §3-5: Piecewise window optimization...")
    alpha_candidates = np.linspace(0.1, 0.9, 30)
    windows = [(-5, -2), (-2, -0.5), (-0.5, 0.5), (0.5, 2), (2, 5)]

    for stype, U, V in [('quadratic', 3.0, 0), ('cubic', 3.0, 1.5), ('multipole', 2.5, 0.8)]:
        window_args = []
        for wlo, whi in windows:
            mask = (omega_grid_200 >= wlo) & (omega_grid_200 < whi)
            omega_slice = omega_grid_200[mask]
            window_args.append((omega_slice, eps0, stype, U, V, alpha_candidates))

        futures_win = client.map(_piecewise_window_sweep, window_args)
        win_alphas = client.gather(futures_win)
        np.savez(os.path.join(CACHE_DIR, f'piecewise_{stype}.npz'),
                 window_alphas=np.array(win_alphas),
                 windows=np.array(windows))
        print(f"    -> Saved piecewise_{stype}.npz")

    # ============================================================
    # 3. RP RESONANCE NOTEBOOK
    # ============================================================
    print("\n" + "=" * 60)
    print("3. RP resonance notebook computations")
    print("=" * 60)

    # §3: Ulam baseline
    print("  §3: Ulam baseline (N=50..1600)...")
    N_ulam_values = [50, 100, 200, 400, 800, 1200, 1600]
    futures_ulam = client.map(_ulam_eigenvalues, N_ulam_values)
    ulam_results = client.gather(futures_ulam)
    ulam_dict = {N: eigs for N, eigs in zip(N_ulam_values, ulam_results)}

    np.savez(os.path.join(CACHE_DIR, 'ulam_baseline.npz'),
             N_values=np.array(N_ulam_values),
             eigenvalues=np.array(ulam_results))
    print("  -> Saved ulam_baseline.npz")

    # §4: Weighted Chebyshev at various N
    print("  §4: Weighted Chebyshev (N=8..96)...")
    N_cheb_values = [8, 16, 24, 32, 48, 64, 96]

    def _cheb_eigs(N):
        L = _chebyshev_transfer_weighted(N, max(500, 5 * N))
        eigs = np.linalg.eigvals(L)
        return np.sort(np.abs(eigs))[::-1][:6]

    futures_cheb = client.map(_cheb_eigs, N_cheb_values)
    cheb_results = client.gather(futures_cheb)

    np.savez(os.path.join(CACHE_DIR, 'cheb_transfer.npz'),
             N_values=np.array(N_cheb_values),
             eigenvalues=np.array(cheb_results))
    print("  -> Saved cheb_transfer.npz")

    # §6: Other parameter values
    print("  §6: Other parameter values...")
    a_values_rp = [3.57, 3.7, 3.9, 4.0]
    N_cheb_test = 48
    N_ulam_test = 800

    rp_other = {}
    for a_val in a_values_rp:
        L_u, _ = frobenius_perron_matrix(logistic, a_val, N=N_ulam_test)
        eigs_u = np.sort(np.abs(np.linalg.eigvals(L_u)))[::-1][:10]
        try:
            L_c = _chebyshev_transfer_weighted(N_cheb_test, 1000, a_param=a_val)
            eigs_c = np.sort(np.abs(np.linalg.eigvals(L_c)))[::-1][:10]
        except Exception:
            eigs_c = np.zeros(10)
        rp_other[a_val] = (eigs_u, eigs_c)

    np.savez(os.path.join(CACHE_DIR, 'rp_other_params.npz'),
             a_values=np.array(a_values_rp),
             ulam_eigs=np.array([rp_other[a][0] for a in a_values_rp]),
             cheb_eigs=np.array([rp_other[a][1] for a in a_values_rp]))
    print("  -> Saved rp_other_params.npz")

    # §7: Jacobi basis comparison
    print("  §7: Jacobi basis comparison...")
    bases = {
        'Chebyshev T (arcsine)': (-0.5, -0.5),
        'Legendre (uniform)': (0.0, 0.0),
        'Chebyshev U': (0.5, 0.5),
        'Jacobi (1, 1)': (1.0, 1.0),
        'Jacobi (-0.25, -0.25)': (-0.25, -0.25),
    }
    N_test_j = 48
    jacobi_args = [(N_test_j, aj, bj, 800, 4.0) for aj, bj in bases.values()]
    futures_jacobi = client.map(_jacobi_transfer_weighted, jacobi_args)
    jacobi_mats = client.gather(futures_jacobi)

    basis_names = list(bases.keys())
    basis_eigs = []
    for L_j in jacobi_mats:
        eigs = np.sort(np.abs(np.linalg.eigvals(L_j)))[::-1][:6]
        basis_eigs.append(eigs)

    np.savez(os.path.join(CACHE_DIR, 'jacobi_comparison.npz'),
             basis_names=np.array(basis_names),
             eigenvalues=np.array(basis_eigs))
    print("  -> Saved jacobi_comparison.npz")

    # §8: Convergence rate study
    print("  §8: Convergence rate study...")
    N_conv = [8, 12, 16, 20, 24, 32, 40, 48, 56, 64, 80, 96]
    conv_args = [(N, -0.5, -0.5, max(500, 5 * N), 4.0) for N in N_conv]
    futures_conv = client.map(_jacobi_transfer_weighted, conv_args)
    conv_mats = client.gather(futures_conv)

    conv_eigs = []
    for L_c in conv_mats:
        eigs = np.sort(np.abs(np.linalg.eigvals(L_c)))[::-1][:4]
        conv_eigs.append(eigs)

    np.savez(os.path.join(CACHE_DIR, 'convergence_rate.npz'),
             N_values=np.array(N_conv),
             eigenvalues=np.array(conv_eigs))
    print("  -> Saved convergence_rate.npz")

    # ============================================================
    # Done
    # ============================================================
    elapsed = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"All precomputation complete in {elapsed:.1f}s")
    print(f"Results in: {CACHE_DIR}")
    print(f"{'=' * 60}")

    client.close()
    cluster.close()


if __name__ == '__main__':
    precompute_all()
