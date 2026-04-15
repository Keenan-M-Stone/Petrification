# Petrification

**Fixed-point dynamics applied to eigenvalue problems, chaos control, and perturbation detection.**

Author: Keenan M Stone  
Origin: May 2017 (original scripts) — April 2026 (current notebook suite)

---

## What This Is

Petrification explores what we'll call the **alpha-transform**, proposed as a simple and obvious way to stabilize, destabilize, and reverse the stability of fixed points in cobweb diagrams — the map

$$g(x) = \alpha f(x) + (1-\alpha)x$$

which preserves the fixed points of $f$ while rescaling their stability. The project investigates what this simple operation reveals when applied systematically to dynamical systems, quantum mechanics, and many-body physics.

The transform is part of the Krasnoselskii-Mann family of relaxation methods. What appears to be new — based on a literature search conducted April 2026 (see [alpha_transform.ipynb §10.5](notebooks/alpha_transform.ipynb) for details) — is the systematic application of position-dependent $\alpha(x)$ as a diagnostic and control tool. No prior work was found using variable relaxation profiles to detect or reconstruct unknown perturbations.

## Key Results

### Established (validated, not novel)

- **Stability inversion theorems** (Props 1-5): clean characterization of how $\alpha$ rescales, inverts, and zeroes fixed-point derivatives
- **$\alpha(x)$ generalization**: position-dependent relaxation with independent control at each fixed point; the $\alpha'$ term vanishes at fixed points
- **Projective eigenstate interpretation**: eigenstates of $\hat{A}$ are fixed points of the projective action $\Pi_A$ on $\mathbb{P}(\mathcal{H})$
- **Frobenius-Perron spectrum**: transfer operator eigenvalues of the logistic map reproduce known Koopman results

### Dead ends (tested, found lacking)

- **Riccati-Bloch eigenvalue computation**: 1200-45000× slower than matrix diagonalization
- **Alpha-transform for eigenstate selection**: trades sensitivity for stability — wrong tradeoff
- **Effective potential from invariant measure**: boundary singularities break self-consistency

### Novel findings

| Finding | Where | Significance |
|---------|-------|-------------|
| **$\alpha^* \approx 0.5$ universality for scalar Dyson equations** | [crossover_alpha_turbiner](notebooks/crossover_alpha_turbiner.ipynb) | Constant mixing $\lambda = 0.5$ is used empirically in GW codes (Chibani et al. 2016) but never derived. Our analytical derivation from fixed-point stability is novel. |
| **Perturbation detection via $\alpha(x)$ profiles** | [perturbation_detection](notebooks/perturbation_detection.ipynb), [quantum](notebooks/quantum_perturbation_detection.ipynb) | $V'_\text{pert} = kx(\alpha(x) - 1)$ exactly reconstructs unmodeled forces. No prior art found (April 2026 search). |
| **Static chaos control via $\alpha(x)$** | [inverse_correspondence](notebooks/inverse_correspondence.ipynb) | Localized $\alpha(x)$ converts chaotic logistic map to superstable fixed point. Related to OGY but as static map modification. |
| **Lyapunov reflection symmetry $\Lambda(\alpha) = \Lambda(-\alpha)$** | [lyapunov_alpha_relationship](notebooks/lyapunov_alpha_relationship.ipynb), [inverse_correspondence](notebooks/inverse_correspondence.ipynb) | $g_\alpha$ and $g_{-\alpha}$ are metrically conjugate ($\mathbb{Z}_2$ symmetry of the Krasnoselskii-Mann family). Reduces $\|\alpha\|$ genuinely stabilizes: $\Lambda(0.5) = -0.74$ at $a=4$ vs $\Lambda(1) = +0.69$. |
| **α(x) dominates OGY/Pyragas for chaos control** | [chaos_control_comparison](notebooks/chaos_control_comparison.ipynb) | Basin of attraction 200/200 vs OGY 43/200. Static map modification is fundamentally more robust than feedback methods. |
| **Z₂ symmetry breaks for complex Dyson iteration** | [kadanoff_baym_alpha](notebooks/kadanoff_baym_alpha.ipynb) | $\Lambda(\alpha) = \Lambda(-\alpha)$ is specific to real-valued maps. Negative α completely fails for Kadanoff-Baym/Dyson (0/200 converged). |
| **α-relaxation improves energy conservation in KB equations** | [kadanoff_baym_alpha](notebooks/kadanoff_baym_alpha.ipynb) | Energy drift decreases monotonically with α: 0.195 → 0.091 → 0.058 for α = 1.0, 0.5, 0.3. α-relaxation regularizes the 2nd Born self-energy. |

## Notebooks

Listed in recommended reading order.

### 1. [alpha_transform.ipynb](notebooks/alpha_transform.ipynb) — Core Theory

The foundation notebook. Develops the alpha-transform from visual intuition (subtract bisector → scale → add back) through formal proofs (Props 1-5) to the position-dependent generalization $\alpha(x)$. Includes quantum eigenstate demonstrations and literature placement.

**Sections:** Literature review → Visual discovery → Formal theory → $\alpha(x)$ generalization → Literature comparison → Demonstrations → Assessment

### 2. [eigenstate_fixedpoint_correspondence.ipynb](notebooks/eigenstate_fixedpoint_correspondence.ipynb) — Full Investigation

The largest notebook (~3300 lines). Systematically tests the eigenstate-fixed-point correspondence with prediction-first methodology. Includes 11 experiments across five directions (forward QM, alpha for eigenstates, reverse spectral→chaos, conceptual quantization, RP resonances).

**Key experiments:** Riccati-shooting benchmark, alpha-relaxed Riccati, Frobenius-Perron operator, effective potential, Ruelle-Pollicott resonances, Riccati transfer operator

**Honest verdict:** Most directions are dead ends or known territory. Value lies in conceptual insight and in systematically closing off unproductive avenues.

### 3. [turbiner_riccati_bloch.ipynb](notebooks/turbiner_riccati_bloch.ipynb) — Riccati-Bloch Deep Dive

Focused treatment of Turbiner's nonlinearization: the Riccati-Bloch equation $y' = y^2 - 2V + 2E$ has bounded trajectories if and only if $E$ is an eigenvalue. Detailed trajectory analysis and boundary condition experiments.

### 4. [crossover_alpha_turbiner.ipynb](notebooks/crossover_alpha_turbiner.ipynb) — Dyson Equation (Novel)

Applies alpha-relaxation to the Dyson equation $G = 1/(\omega - \varepsilon_0 - \Sigma[G])$. Discovers that constant $\alpha \approx 0.5$ provides universal convergence across all coupling strengths. Diagnoses the bootstrap failure of adaptive $\alpha(\omega)$. Connects to the SCN/Nullified framework. Includes matrix Dyson and DMFT future directions.

**This is the strongest candidate for publication** (Phys. Rev. B or Comp. Phys. Comm.).

### 5. [inverse_correspondence.ipynb](notebooks/inverse_correspondence.ipynb) — Regular vs. Inverted Dynamics

Exploratory notebook comparing regular ($\alpha = 1$) and inverted ($\alpha = -1$) logistic map dynamics. Dual bifurcation diagrams, period analysis, Lyapunov comparison, continuous $\alpha$ sweep, full $(\alpha, a)$ phase diagram, and position-dependent $\alpha(x)$ experiments. Includes connection to OGY chaos control.

### 6. [perturbation_detection.ipynb](notebooks/perturbation_detection.ipynb) — Classical Perturbation Detection

Demonstrates that measuring the $\alpha(x)$ profile of a damped oscillator exactly reconstructs perturbation forces: parallel Hooke springs, displaced equilibria, quartic anharmonicity. Blind reconstruction test.

### 7. [quantum_perturbation_detection.ipynb](notebooks/quantum_perturbation_detection.ipynb) — Quantum Perturbation Detection

Extends perturbation detection to the hydrogen atom. Sequential stacking procedure: Coulomb baseline → detect fine structure (1/r³) → detect vacuum polarization (Yukawa). Shape of $\alpha(r) - 1$ discriminates power-law vs. exponential perturbations. Includes blind Gaussian well test.

### Short-term Investigations

These notebooks test specific predictions from the core theory with Dask-parallelized precomputation (cache in `notebooks/cache/`).

#### 8. [matrix_dyson_universality.ipynb](notebooks/matrix_dyson_universality.ipynb) — Matrix Dyson $\alpha^*$

Tests whether $\alpha^* \approx 0.5$ holds for multi-orbital Dyson equations. Result: holds robustly for diagonal and hybridized self-energies (stays in [0.43, 0.55]); drifts to ~0.62–0.65 for full matrix Σ at higher dimensions, but α=0.5 still achieves 100% convergence.

#### 9. [lyapunov_alpha_relationship.ipynb](notebooks/lyapunov_alpha_relationship.ipynb) — Lyapunov Reflection Symmetry

Discovers the $\mathbb{Z}_2$ reflection $\Lambda(\alpha) = \Lambda(-\alpha)$ — regular and inverted maps have identical Lyapunov exponents ($r = 0.9999$ across 500 parameter values). Reducing $|\alpha|$ below 1 genuinely stabilizes: $\Lambda(0.5) = -0.74$ at $a=4$ vs $\Lambda(1) = +0.69$.

#### 10. [piecewise_alpha_dyson.ipynb](notebooks/piecewise_alpha_dyson.ipynb) — Piecewise $\alpha(\omega)$ (Negative)

Tests whether frequency-dependent $\alpha$ improves Dyson convergence. Result: $\alpha^*(\omega)$ barely varies (range [0.503, 0.528]), and piecewise optimization just rediscovers $\alpha \approx 0.5$ everywhere.

#### 11. [rp_resonance_basis.ipynb](notebooks/rp_resonance_basis.ipynb) — Ruelle-Pollicott Convergence (Negative)

Tests Jacobi polynomial bases for RP resonance computation. Neither Ulam nor weighted Chebyshev converges to exact sub-leading eigenvalues. Convergence rate is $O(N^{0.01})$ — essentially flat.

### Medium-term Investigations

These notebooks require new frameworks beyond the core theory.

#### 12. [chaos_control_comparison.ipynb](notebooks/chaos_control_comparison.ipynb) — α(x) vs OGY vs Pyragas

Head-to-head comparison of three chaos control methods. α(x) dominates: 200/200 basin of attraction vs OGY 43/200, ~2× more noise-robust than Pyragas, and faster convergence. The key insight: α-control is a structural map modification (no feedback loop), making it fundamentally more robust.

#### 13. [experimental_perturbation_detection.ipynb](notebooks/experimental_perturbation_detection.ipynb) — Experimental Detection Limits

Tests α(x) perturbation reconstruction under realistic experimental conditions: optical trap simulation, noise robustness (works down to SNR ≈ 3 dB), multi-perturbation stacking, perturbation classification from α-profile shape, and Shannon channel capacity analysis.

#### 14. [kadanoff_baym_alpha.ipynb](notebooks/kadanoff_baym_alpha.ipynb) — Kadanoff-Baym / Non-equilibrium Green's Functions

Applies α-relaxation to the two-time Dyson equation on the Keldysh contour. Three key findings: (1) optimal $(\alpha_R, \alpha_<) = (0.10, 0.10)$ — lower than scalar Dyson's $\alpha^* \approx 0.5$; (2) Z₂ symmetry breaks for complex-valued Dyson iteration; (3) smaller α improves energy conservation (drift 0.058 at α=0.3 vs 0.195 at α=1).

### Cross-Theory Connections

[cross_theory_connections.md](notebooks/cross_theory_connections.md) catalogs connections to other mathematical frameworks (group theory, gauge theory, information theory, graph theory, tensor networks, fluid mechanics, etc.) that may simplify, generalize, or provide alternative notation for the alpha-transform results. Priority items: deviation-space coordinates (§12.2), Jacobi polynomial basis for improved RP resonances (§10.1), and information-theoretic bounds on perturbation detection (§3.2).

## Python Package: `petrification/`

```
petrification/
├── __init__.py          # Package docstring and module listing
├── maps.py              # Discrete maps: logistic, Migdal-Kadanoff, exponential saturation
├── transforms.py        # Alpha-transform, optimal alpha, alpha(x) construction
├── iteration.py         # iterate(), iterate_transformed(), cobweb_data(), find_fixed_points()
├── bifurcation.py       # compute_bifurcation(), compute_bifurcation_transformed()
├── potentials.py        # Quantum potentials: harmonic, anharmonic, double-well, Morse, Coulomb
├── quantum.py           # Hamiltonian discretization, Riccati-Bloch, spectral scanning,
│                        #   power iteration, Numerov, Frobenius-Perron
├── dyson.py             # Scalar/matrix/lattice Dyson equation iteration
└── oscillators.py       # Damped oscillator simulation, fitting, perturbation reconstruction
```

### Key Functions

| Module | Function | Purpose |
|--------|----------|---------|
| `maps` | `logistic(a, x)` | Standard logistic map $f(x) = ax(1-x)$ |
| `transforms` | `alpha_transform(f, alpha, a, x)` | Apply $g = \alpha f + (1-\alpha)x$ |
| `transforms` | `compute_optimal_alpha(f, a, x_domain, boundary)` | Find $\alpha^* = 1/(1 - f'_{\max})$ |
| `transforms` | `make_alpha_func(f', a, smooth, cap)` | Build callable $\alpha(x)$ with singularity capping |
| `iteration` | `iterate_transformed(f, alpha, a, x0, n_iter)` | Iterate with constant or callable $\alpha$ |
| `iteration` | `cobweb_data(f, a, x0, n_iter, alpha)` | Generate cobweb diagram coordinates |
| `bifurcation` | `compute_bifurcation(f, a_range, ...)` | Bifurcation diagram: long-term attractor points |
| `bifurcation` | `compute_bifurcation_transformed(f, alpha, ...)` | Bifurcation under alpha-transform |
| `quantum` | `solve_eigenstates(V, x_grid, n_states)` | Diagonalize discretized Hamiltonian |
| `quantum` | `riccati_solve(V, E, x_span, y0)` | Integrate Riccati-Bloch equation |
| `quantum` | `frobenius_perron_matrix(f, a, N, x_range)` | Transfer operator via Monte Carlo/Ulam |
| `dyson` | `scalar_dyson_iterate(...)` | Alpha-relaxed Dyson equation iteration |
| `oscillators` | `simulate(V', gamma, m, x0, v0, ...)` | RK45 damped oscillator integration |
| `oscillators` | `infer_perturbation(...)` | Reconstruct $V'_\text{pert}$ from $\alpha(x)$ profile |

## Standalone Scripts (`scripts/`)

| Script | Purpose |
|--------|--------|
| `scripts/Bifurcation_Transform.py` | Bifurcation diagrams for alpha-transformed logistic map |
| `scripts/Bifurcation_logistic_gif.py` | Animated GIF of bifurcation across $\alpha$ values |
| `scripts/Iteration_Transform.py` | Cobweb diagram visualization with alpha-transform |
| `scripts/eigen_fixedpoint.py` | Eigenstate-fixed-point comparison (original 2017 script) |
| `scripts/eigen_fixedpoint_v2.py` | Enhanced eigenstate analysis with visualization |

## Open Questions and Future Directions

### Short-term (directly actionable)

- [x] **Matrix Dyson equations**: Does $\alpha^* \approx 0.5$ hold for multi-orbital Dyson in DMFT? This is the natural next step from the scalar result.
    - **Has substance**: α* ≈ 0.5 holds robustly for diagonal and hybridized self-energies (stays in [0.43, 0.55]). Full matrix Σ drifts to ~0.62-0.65 at higher dimensions, but α=0.5 still achieves 100% convergence (200/200 frequencies vs naive 42/200). The superoperator analysis gives α*=0.520 for d=3. There's a real story here about when and why the universality breaks.
- [x] **Lyapunov exponent relationship**: How do Lyapunov exponents of regular and $\alpha$-transformed maps relate? Run §4 of inverse_correspondence and characterize the $(\alpha, a)$ phase diagram.
    - **Has substance (with surprise)**: The Lyapunov exponent satisfies a **reflection symmetry** $\Lambda(\alpha) = \Lambda(-\alpha)$ — confirmed at individual $\alpha$ values ($\Lambda(+0.5) = \Lambda(-0.5)$ exactly) and in the $(\alpha, a)$ phase diagram (chaos boundaries symmetric around $\alpha = 0$). The $\alpha = \pm 1$ case is the most physically meaningful: regular vs inverted maps have $r = 0.9999$ across 500 parameter values, $\text{RMSE} = 0.003$. However, $\Lambda(\alpha) \neq \Lambda(1)$ for general $|\alpha| \neq 1$ — at $a=4$, $\Lambda(0.5) = -0.74$ while $\Lambda(1) = +0.69$, so reducing $|\alpha|$ below 1 genuinely stabilizes the dynamics. The symmetry is $\mathbb{Z}_2$ (sign flip), not full invariance under the group action.
- [x] **Piecewise $\alpha(\omega)$**: Can we overcome the bootstrap failure of adaptive $\alpha$ by using different constant $\alpha$ values in different frequency windows?
  - **Negative**: α*(ω) barely varies across frequency: quadratic [0.504, 0.525], cubic [0.503, 0.528], multipole [0.535, 0.542]. Piecewise optimization just rediscovers α ≈ 0.5 everywhere. It's only better than a bad constant choice (piecewise saves 57 iters vs α=0.9, but α=0.55 already gets 55.3 avg iters). The window count sweep confirms: increasing from 1 to 20 windows barely changes anything (55.3 → 55.1 avg iters). Not worth pursuing — a single global α ≈ 0.5 is already near-optimal for scalar Dyson.
- [x] **Minimal discretization for RP resonances**: What resolution is needed for exact sub-leading Ruelle-Pollicott eigenvalues? Consider Jacobi polynomial basis weighted by arcsine measure (see [cross-theory connections §10](notebooks/cross_theory_connections.md)).
    - **Negative**: Neither Ulam nor weighted Chebyshev converges to the exact RP eigenvalues. Ulam at N=1600 still gives |λ₁| ≈ 0.57 instead of 0.50. Chebyshev is worse — |λ₀| drifts to 1.79 and |λ₁| error plateaus at 0.48 (essentially no convergence). The Jacobi basis comparison shows Legendre performs best for |λ₀| (0.9995 vs exact 1.0) but all bases fail badly on |λ₁|. Convergence rate is O(N^0.01) — effectively flat. The chebyshev_transfer_weighted implementation likely has a bug in the quadrature or weight handling, and even if fixed, the Bernoulli map's singular invariant measure makes spectral methods struggle.

### Medium-term (requires new framework)

- [x] **$\alpha(x)$ perturbation detection in experiment**: Apply the reconstruction formula $V'_\text{pert} = kx(\alpha(x) - 1)$ to real experimental data (e.g., optical trap perturbation measurement)
    - **Has substance**: Noise robustness confirmed — α(x) reconstruction works down to SNR ≈ 3 dB with useful fidelity. Optical trap simulation demonstrates full perturbation detection pipeline (parallel Hooke spring, displaced equilibrium, quartic anharmonicity). Perturbation classification from α-profile shape achieved with 85%+ accuracy. In a 15×15 amplitude/noise sweep, detection threshold follows a clean contour in (amplitude, noise) space. Multi-perturbation stacking works: sequential α-profile subtraction isolates individual perturbation signatures. Shannon channel capacity analysis quantifies the information content of α(x) as a diagnostic signal.
- [x] **Chaos control comparison**: Quantitative comparison of $\alpha(x)$ chaos suppression vs. OGY and Pyragas methods — compute basins of attraction, robustness to noise, control effort
    - **Has substance**: α(x) dominates OGY and Pyragas across all metrics. Basins of attraction: α(x) stabilizes from 200/200 initial conditions at $a=3.8$ vs OGY 43/200 and Pyragas 112/200. Noise robustness: α(x) maintains fixed-point convergence at noise levels where OGY fails (~2× more robust). Convergence speed: α(x) reaches steady state in ~20 iterations vs 50+ for Pyragas. Control effort: α(x) pays for this by reshaping the entire map (higher total effort), but per-iterate effort is lower because no real-time orbit detection is needed. The Lyapunov reflection symmetry $\Lambda(\alpha) = \Lambda(-\alpha)$ is confirmed: $\alpha = -0.5$ and $\alpha = +0.5$ give identical stabilization. Period-2 orbit stabilization also demonstrated. The key insight: α-control is a structural modification (no feedback loop), making it fundamentally more robust than OGY/Pyragas feedback methods.
- [x] **Non-equilibrium Green's functions**: Alpha-relaxation for the Kadanoff-Baym equations (two-time Dyson equation on the Keldysh contour)
    - **Has substance (with caveats)**: Three key results: (1) α-relaxation accelerates KB self-consistency convergence — optimal at $(\alpha_R, \alpha_<) = (0.10, 0.10)$ with residual $3.88 \times 10^{-3}$ vs much larger at naive $\alpha = 1$, confirming that the scalar Dyson $\alpha^* \approx 0.5$ does **not** transfer directly (KB optimal is lower). (2) **Z₂ symmetry breaks**: $\alpha = -0.5$ completely fails to converge (0/200 frequencies) while $\alpha = +0.5$ converges perfectly — the Lyapunov reflection symmetry $\Lambda(\alpha) = \Lambda(-\alpha)$ is specific to real-valued maps and does not carry over to complex-valued Dyson iteration. (3) **α-relaxation improves energy conservation**: energy drift decreases monotonically with decreasing α (drift = 0.195 at α=1.0, 0.091 at α=0.5, 0.058 at α=0.3), showing that α-relaxation acts as a regularizer preserving the conserving properties of the 2nd Born self-energy approximation. Caveat: the KB solver uses explicit time-stepping with self-energy under-relaxation (η=0.3) and amplitude clamping, which limits the U range and time extent accessible. The qualitative conclusions are robust but quantitative optima will shift with more sophisticated solvers.

### Speculative (may or may not lead anywhere)

- [ ] **Gauge-theoretic formulation**: Does treating $\alpha(x)$ as a gauge field (with the operator group $\Gamma_\alpha \Gamma_\beta = \Gamma_{\alpha\beta}$) produce conservation laws or topological invariants?
    - **Prediction (partially supported, narrow scope)**: The Lyapunov invariance $\Lambda(\alpha) = \Lambda(-\alpha)$ is a conservation law of the $\Gamma_\alpha$ group under the $\mathbb{Z}_2$ reflection $\alpha \mapsto -\alpha$. The piecewise Dyson result ($\alpha^*$ nearly constant across $\omega$) means the gauge field is "flat" in the Dyson context — zero curvature, no interesting gauge physics. However, the matrix Dyson result ($\alpha^*$ drifts from 0.5 to 0.62-0.65 as orbital dimension increases) suggests non-trivial gauge curvature in orbital space. Prediction: the gauge formulation is unproductive for scalar problems but becomes interesting for multi-orbital systems where the optimal mixing matrix $\boldsymbol{\alpha}$ has off-diagonal structure.
- [ ] **Concavity inversion in optimization**: Since $g'' = \alpha f''$, negative $\alpha$ flips convexity. Can this be used to escape local minima in non-convex optimization by temporarily inverting the landscape?
    - **Prediction (unlikely to be useful)**: The Lyapunov data shows $\Lambda(\alpha) = \Lambda(-\alpha)$, so negative $\alpha$ gives exactly the same global convergence rate as positive $\alpha$ of the same magnitude. Concavity flipping doesn't escape local minima any better than standard relaxation — it just accesses the same dynamics through a different route. The period diagram does show different local attractor structures at $\alpha$ vs $-\alpha$ (different periodic orbits become stable), so in principle different fixed points are reachable. But this is redundant with simply varying $|\alpha|$ to change basin boundaries, which is standard relaxation parameter tuning.
- [ ] **Koopman operator connection**: Does the $\alpha$-transform correspond to a similarity transformation on the Koopman operator? If so, the spectral invariance (equal Lyapunov exponents) follows immediately.
    - **Prediction (partially true — $\mathbb{Z}_2$ not full group)**: The Lyapunov data reveals the situation is more specific than originally conjectured. Lyapunov exponents satisfy $\Lambda(\alpha) = \Lambda(-\alpha)$ exactly (the reflection symmetry), and $\Lambda(1) = \Lambda(-1)$ to high precision ($r = 0.9999$ across 500 parameter values). But $\Lambda(\alpha) \neq \Lambda(1)$ for general $|\alpha| \neq 1$ — at $a=4$, $\Lambda(0.5) = -0.74$ while $\Lambda(1) = +0.69$. So the Koopman similarity holds only under $\alpha \mapsto -\alpha$ (a $\mathbb{Z}_2$ symmetry), not the full multiplicative group. The Koopman operators of $g_\alpha$ and $g_{-\alpha}$ should be isospectral, but $g_{0.5}$ and $g_1$ have genuinely different spectra. This is provable for maps with the involution symmetry $f(1-x) = f(x)$ (like the logistic map), since $g_{-\alpha}(1-x) = 1 - g_\alpha(x)$ in that case. The clean $\mathbb{Z}_2$ result is likely a short paper on its own.
- [ ] **Renormalization group connection**: The Migdal-Kadanoff map is already an RG transformation. What does $\alpha$-transforming an RG map mean physically?
    - **Prediction (conceptually clear, probably not new)**: Alpha-transforming an RG map with $|\alpha| < 1$ slows the RG flow without changing its fixed points — this is exactly what alpha-relaxation does in the Dyson/DMFT context (the DMFT self-consistency loop IS an RG-like iteration). The $\alpha^* \approx 0.5$ result translates to: the optimal RG flow speed is about half the naive rate. This is real but may not lead to new physics beyond what the Dyson results already show. It becomes interesting if different universality classes require different $\alpha^*$, which would mean the optimal relaxation rate encodes the RG fixed-point structure.

## Environment

```bash
conda activate moonstone
# Python 3.13+, numpy, scipy, matplotlib, sympy, dask
```

## Git History (highlights)

```
7ac7a7d  Add inverse correspondence notebook
4de2484  Add quantum perturbation detection notebook
b6f6d5f  Add perturbation detection experiment
ba87eba  Add motivational sections to alpha_transform notebook
519eb31  Split monolithic notebook into three focused notebooks
45d8b6c  Add Dyson equation α-relaxation experiments
64a8b46  Add Coulomb potential, Chebyshev/Ulam transfer operators
e5db868  α(x) generalization: position-dependent transform
41ea358  Exp 4-6: Dask-parallelized RP resonances
59e18ad  Prediction-first experiments: Riccati, Frobenius-Perron
42659fd  Initial commit
```

## Afterword

If any of this ever finds any use to anyone, then I'd like to go ahead and acknowledge my
grad school family that continued to support me long after the nightmare ended.
With special thanks to:

- Amara Katabarwa, for his thoughtful insights into the original project.
- Eric Suter, for his encouragements and for coining the term "Petrification".
- Brandon Campbell, who was also in that class with us.
- Antonio Mantica, who may have been in that class with us.
- Lauren Sgro, did you take that class with us?

Dedicated to the memory of Howard Lee.  
"This began with a piece of chalk and heaven opened up." - H. Lee
