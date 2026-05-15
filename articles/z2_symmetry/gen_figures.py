"""Generate figures for z2_symmetry paper."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Vectorized Lyapunov computation ──────────────────────────────────────────
def lyap_vec(theta, a, x0=0.4, n=5000, nd=500):
    """Compute Lyapunov exponent for arrays of (theta, a) pairs."""
    # theta, a can be scalars or arrays (broadcast together)
    th = np.asarray(theta, dtype=float)
    av = np.asarray(a,     dtype=float)
    x = np.full_like(th, x0)
    alive = np.ones_like(th, dtype=bool)

    for _ in range(nd):
        fx = av * x * (1 - x)
        x = th * fx + (1 - th) * x
        alive &= (x > 0) & (x < 1)

    s = np.zeros_like(th)
    for _ in range(n):
        fp = av * (1 - 2*x)
        gp = 1 + th * (fp - 1)
        s += np.where(alive & (np.abs(gp) > 1e-30),
                      np.log(np.maximum(np.abs(gp), 1e-30)), 0.0)
        fx = av * x * (1 - x)
        x = th * fx + (1 - th) * x
        alive &= (x > 0) & (x < 1)

    result = np.where(alive, s / n, np.nan)
    return result


# ── Figure 1: Lyapunov vs theta, four panels ─────────────────────────────────
theta_range = np.linspace(-2, 2, 401)
a_test  = [3.2, 3.57, 3.83, 4.0]
labels  = ['$a=3.2$ (period-2)', '$a=3.57$ (chaos onset)',
           '$a=3.83$ (period-3)', '$a=4.0$ (fully chaotic)']
colors  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

fig, axes = plt.subplots(2, 2, figsize=(3.375*2, 2.6*2), constrained_layout=True)
axes = axes.flatten()

for idx, (a, lab, col) in enumerate(zip(a_test, labels, colors)):
    lv = lyap_vec(theta_range, np.full_like(theta_range, a))
    ax = axes[idx]
    ax.plot(theta_range, lv, color=col, lw=1.0)
    ax.axhline(0,  color='k',    lw=0.6)
    ax.axvline( 1, color='gray', ls=':', lw=0.8)
    ax.axvline(-1, color='gray', ls=':', lw=0.8)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-4, 2)
    ax.set_xlabel(r'$\theta$', fontsize=9)
    ax.set_ylabel(r'$\Lambda(\theta)$', fontsize=9)
    ax.set_title(lab, fontsize=9)
    ax.tick_params(labelsize=8)
    panel = chr(ord('a') + idx)
    ax.text(0.04, 0.95, f'({panel})', transform=ax.transAxes,
            fontsize=9, va='top', fontweight='bold')

fig.savefig('figures/lyap_vs_theta.pdf', bbox_inches='tight')
print('lyap_vs_theta.pdf saved')

# ── Figure 2: (theta, a) phase diagram — vectorized row-by-row ───────────────
NTH = 120
NA  = 120
theta_2d = np.linspace(-2, 2, NTH)
a_2d     = np.linspace(2.5, 4.0, NA)
lyap_map = np.full((NTH, NA), np.nan)

for i, th in enumerate(theta_2d):
    lyap_map[i, :] = lyap_vec(np.full(NA, th), a_2d, n=3000, nd=300)

fig2, ax2 = plt.subplots(figsize=(3.375, 2.6))
extent = [a_2d[0], a_2d[-1], theta_2d[0], theta_2d[-1]]
lm_clip = np.clip(lyap_map, -3, 1.5)
im = ax2.imshow(lm_clip, extent=extent, aspect='auto', origin='lower',
                cmap='RdBu_r', vmin=-3, vmax=1.5)
plt.colorbar(im, ax=ax2, label=r'$\Lambda(\theta, a)$',
             fraction=0.046, pad=0.04)
lm_filled = np.where(np.isnan(lyap_map), 10, lyap_map)
ax2.contour(a_2d, theta_2d, lm_filled, levels=[0],
            colors='black', linewidths=1.0)
ax2.axhline( 1, color='white', ls='--', lw=0.9, label=r'$\theta=+1$')
ax2.axhline(-1, color='white', ls=':', lw=0.9, label=r'$\theta=-1$')
ax2.set_xlabel('Parameter $a$', fontsize=9)
ax2.set_ylabel(r'$\theta$', fontsize=9)
ax2.legend(fontsize=7, frameon=False, loc='upper left')
ax2.tick_params(labelsize=8)
fig2.tight_layout(pad=0.4)
fig2.savefig('figures/phase_diagram.pdf', bbox_inches='tight')
print('phase_diagram.pdf saved')
