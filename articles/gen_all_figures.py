"""Generate all figures for perturbation_probing and chaos_control papers."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path('../../..').resolve()))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ── logistic helpers ──────────────────────────────────────────────────────────
def logistic(a, x):
    return a * x * (1.0 - x)

def logistic_fp(a):
    return 1.0 - 1.0 / a

def logistic_deriv(a, x):
    return a * (1.0 - 2.0 * x)

# ─────────────────────────────────────────────────────────────────────────────
#  PERTURBATION PROBING FIGURES
# ─────────────────────────────────────────────────────────────────────────────

A = 3.7
X1 = logistic_fp(A)   # ≈ 0.7297

# ── Fig 1: FP drift vs epsilon (IFT theorem) — two FPs, sine perturbation ────
def f_sine(a, x, eps, omega):
    return a * x * (1.0 - x) + eps * np.sin(omega * x)

def find_fp(a, eps, omega, guess, bracket=(1e-6, 0.9999)):
    g = lambda x: f_sine(a, x, eps, omega) - x
    lo = max(bracket[0], guess - 0.35)
    hi = min(bracket[1], guess + 0.35)
    if g(lo) * g(hi) > 0:
        lo, hi = bracket
    try:
        return brentq(g, lo, hi, xtol=1e-12)
    except ValueError:
        return np.nan

def ift_fp(x_star, eps, omega, a):
    denom = a * (1.0 - 2.0 * x_star) - 1.0
    return x_star - eps * np.sin(omega * x_star) / denom

omega = 10.0
eps_vals = np.linspace(-0.05, 0.05, 201)

fp1_num = np.array([find_fp(A, e, omega, X1) for e in eps_vals])
fp1_ift = np.array([ift_fp(X1, e, omega, A) for e in eps_vals])
fp0_num = np.array([find_fp(A, e, omega, 0.0, (0.0, 0.15)) for e in eps_vals])
fp0_ift = np.array([ift_fp(0.0, e, omega, A) for e in eps_vals])

fig, axes = plt.subplots(1, 2, figsize=(3.375*2, 2.4), constrained_layout=True)
for ax, num, ift, x0, lab in zip(
        axes,
        [fp1_num, fp0_num],
        [fp1_ift, fp0_ift],
        [X1, 0.0],
        [r'Non-trivial FP $x^*_1$', r'Trivial FP $x^*_0$']):
    ax.plot(eps_vals, num, 'b-',  lw=1.2, label='Numerical')
    ax.plot(eps_vals, ift, 'r--', lw=1.2, label='IFT (linear)')
    ax.axhline(x0, color='gray', lw=0.7, ls=':')
    ax.set_xlabel(r'$\varepsilon$', fontsize=9)
    ax.set_ylabel(r'$x^*(\varepsilon)$', fontsize=9)
    ax.set_title(lab, fontsize=9)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)

fig.savefig('perturbation_probing/figures/fp_drift.pdf', bbox_inches='tight')
print('fp_drift.pdf saved')

# ── Fig 2: Null-frequency sweep ───────────────────────────────────────────────
eps_fixed = 0.03
omega_vals = np.linspace(0.5, 60.0, 1200)
denom_1 = A * (1.0 - 2.0 * X1) - 1.0
ift_drift = -(eps_fixed * np.sin(omega_vals * X1)) / denom_1

# spot-check a handful of frequencies
omega_spot = np.array([1.0, 5.0, 10.0, 20.0, 30.0])
fp_spot = np.array([find_fp(A, eps_fixed, w, X1) - X1 for w in omega_spot])

null_k = np.arange(1, int(60 * X1 / np.pi) + 1)
omega_null = null_k * np.pi / X1

fig2, ax2 = plt.subplots(figsize=(3.375, 2.4))
ax2.plot(omega_vals, ift_drift, 'b-', lw=1.0, label=r'IFT: $\Delta x^*_1$')
ax2.scatter(omega_spot, fp_spot, color='red', zorder=5, s=30, label='Numerical')
ax2.axhline(0, color='gray', lw=0.7, ls=':')
nulls_in = omega_null[omega_null < 60]
ax2.vlines(nulls_in, -0.016, 0.016,
           color='green', alpha=0.35, lw=0.7, label=r'Null: $\omega_k = k\pi/x^*_1$')
ax2.set_xlabel(r'$\omega$', fontsize=9)
ax2.set_ylabel(r'$\Delta x^*_1$', fontsize=9)
ax2.set_title(rf'FP drift vs.\ frequency ($\varepsilon={eps_fixed}$, $a={A}$)', fontsize=9)
ax2.legend(fontsize=7, frameon=False)
ax2.tick_params(labelsize=8)
fig2.tight_layout(pad=0.4)
fig2.savefig('perturbation_probing/figures/null_frequencies.pdf', bbox_inches='tight')
print('null_frequencies.pdf saved')

# ── Fig 3: Gaussian selectivity ratio vs sigma ────────────────────────────────
sigma_vals = np.linspace(0.02, 0.20, 200)
X0 = 0.0
ratio = np.exp(-(X1 - X0)**2 / (2 * sigma_vals**2))

fig3, ax3 = plt.subplots(figsize=(3.375, 2.4))
ax3.semilogy(sigma_vals, ratio, 'purple', lw=1.2)
ax3.axvline(0.08, color='gray', ls='--', lw=0.8)
ax3.text(0.082, 1e-15, r'$\sigma=0.08$', fontsize=7, color='gray')
ax3.set_xlabel(r'$\sigma$', fontsize=9)
ax3.set_ylabel(r'$|\Delta x^*_0|/|\Delta x^*_1|$', fontsize=9)
ax3.set_title('Gaussian selectivity ratio (Theorem 3)', fontsize=9)
ax3.tick_params(labelsize=8)
fig3.tight_layout(pad=0.4)
fig3.savefig('perturbation_probing/figures/selectivity.pdf', bbox_inches='tight')
print('selectivity.pdf saved')

# ─────────────────────────────────────────────────────────────────────────────
#  CHAOS CONTROL FIGURES
# ─────────────────────────────────────────────────────────────────────────────

def ogy_iterate(a_nom, x0, n_iter, delta_max=0.05):
    traj = np.zeros(n_iter + 1)
    traj[0] = x0
    x_star = logistic_fp(a_nom)
    fp = logistic_deriv(a_nom, x_star)
    dfda = x_star * (1 - x_star)
    for i in range(n_iter):
        x = traj[i]
        dx = x - x_star
        if abs(dx) < 0.1 and abs(dfda) > 1e-10:
            delta_a = np.clip(-fp * dx / dfda, -delta_max, delta_max)
        else:
            delta_a = 0.0
        traj[i+1] = logistic(a_nom + delta_a, x)
    return traj

def pyragas_opt_K(a):
    fp = logistic_deriv(a, logistic_fp(a))
    K_crit = (1 + fp) / 2
    return np.clip(K_crit - 0.05, -0.95, K_crit - 0.01)

def pyragas_iterate(a, x0, n_iter, tau=1):
    K = pyragas_opt_K(a)
    traj = np.zeros(n_iter + 1)
    traj[0] = x0
    for i in range(min(tau, n_iter)):
        traj[i+1] = logistic(a, traj[i])
    for i in range(tau, n_iter):
        feedback = K * (traj[i - tau] - traj[i])
        traj[i+1] = np.clip(logistic(a, traj[i]) + feedback, 0.001, 0.999)
    return traj

def alpha_ss(a):
    fp = logistic_deriv(a, logistic_fp(a))  # 2 - a
    return 1.0 / (1.0 - fp + 1e-10)        # 1 / (a - 1)

def alpha_func(x, a, sigma=0.1):
    xs = logistic_fp(a)
    ass = alpha_ss(a)
    return 1.0 - (1.0 - ass) * np.exp(-(x - xs)**2 / (2.0 * sigma**2))

def alpha_iterate(a, x0, n_iter, sigma=0.1):
    traj = np.zeros(n_iter + 1)
    traj[0] = x0
    for i in range(n_iter):
        x = traj[i]
        al = alpha_func(x, a, sigma)
        traj[i+1] = al * logistic(a, x) + (1.0 - al) * x
    return traj

# ── Fig 4: Trajectory comparison at a = 3.8 ──────────────────────────────────
a = 3.8
x_star = logistic_fp(a)
n_iter = 100
x0 = 0.2

traj_none = np.zeros(n_iter + 1)
traj_none[0] = x0
for i in range(n_iter):
    traj_none[i+1] = logistic(a, traj_none[i])

traj_ogy = ogy_iterate(a, x0, n_iter)
traj_pyr = pyragas_iterate(a, x0, n_iter)
traj_alp = alpha_iterate(a, x0, n_iter)

fig4, (axT, axE) = plt.subplots(2, 1, figsize=(3.375, 3.6), constrained_layout=True)

axT.plot(range(n_iter+1), traj_none, color='gray', lw=0.5, alpha=0.4, label='Uncontrolled')
axT.plot(range(n_iter+1), traj_ogy, 'b-', lw=1.0, label='OGY')
axT.plot(range(n_iter+1), traj_pyr, 'r-', lw=1.0, label='Pyragas')
axT.plot(range(n_iter+1), traj_alp, 'g-', lw=1.0, label=r'$\alpha(x)$')
axT.axhline(x_star, color='k', ls='--', lw=0.7)
axT.set_ylabel(r'$x_n$', fontsize=9)
axT.set_title(rf'$a={a}$, $x_0={x0}$', fontsize=9)
axT.legend(fontsize=7, frameon=False)
axT.tick_params(labelsize=8)

for traj, lab, col in [(traj_ogy,'OGY','b'), (traj_pyr,'Pyragas','r'), (traj_alp,r'$\alpha(x)$','g')]:
    err = np.abs(traj - x_star)
    err = np.where(err < 1e-14, 1e-14, err)
    axE.semilogy(range(n_iter+1), err, color=col, lw=1.0, label=lab)
axE.set_xlabel('Iteration $n$', fontsize=9)
axE.set_ylabel(r'$|x_n - x^*|$', fontsize=9)
axE.legend(fontsize=7, frameon=False)
axE.tick_params(labelsize=8)

fig4.savefig('chaos_control/figures/trajectories.pdf', bbox_inches='tight')
print('trajectories.pdf saved')

# ── Fig 5: Basin of attraction bar chart ─────────────────────────────────────
a_values = [3.5, 3.7, 3.8, 3.9, 4.0]
x0_grid = np.linspace(0.01, 0.99, 200)
n_basin = 300
tol = 1e-4

basin_ogy = []
basin_pyr = []
basin_alp = []

for a in a_values:
    xs = logistic_fp(a)
    c_o = c_p = c_a = 0
    for x0 in x0_grid:
        if abs(alpha_iterate(a, x0, n_basin)[-1] - xs) < tol:
            c_a += 1
        if abs(ogy_iterate(a, x0, n_basin)[-1] - xs) < tol:
            c_o += 1
        if abs(pyragas_iterate(a, x0, n_basin)[-1] - xs) < tol:
            c_p += 1
    total = len(x0_grid)
    basin_alp.append(100 * c_a / total)
    basin_ogy.append(100 * c_o / total)
    basin_pyr.append(100 * c_p / total)
    print(f'a={a}: α(x)={100*c_a/total:.0f}%, OGY={100*c_o/total:.0f}%, Pyr={100*c_p/total:.0f}%')

x_pos = np.arange(len(a_values))
width = 0.27

fig5, ax5 = plt.subplots(figsize=(3.375, 2.5))
ax5.bar(x_pos - width, basin_alp, width, label=r'$\alpha(x)$', color='g', alpha=0.8)
ax5.bar(x_pos,         basin_ogy, width, label='OGY',          color='b', alpha=0.8)
ax5.bar(x_pos + width, basin_pyr, width, label='Pyragas',      color='r', alpha=0.8)
ax5.set_xticks(x_pos)
ax5.set_xticklabels([f'$a={a}$' for a in a_values], fontsize=8)
ax5.set_ylabel('Basin size (%)', fontsize=9)
ax5.set_ylim(0, 110)
ax5.legend(fontsize=8, frameon=False)
ax5.tick_params(labelsize=8)
fig5.tight_layout(pad=0.4)
fig5.savefig('chaos_control/figures/basins.pdf', bbox_inches='tight')
print('basins.pdf saved')
