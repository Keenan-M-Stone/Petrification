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
| **Lyapunov exponent relationship under constant $\alpha$** | [inverse_correspondence](notebooks/inverse_correspondence.ipynb) | Open question: how do Lyapunov exponents of regular and $\alpha$-transformed maps relate? The $(\alpha, a)$ phase diagram is unexplored. |

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

- [ ] **Matrix Dyson equations**: Does $\alpha^* \approx 0.5$ hold for multi-orbital Dyson in DMFT? This is the natural next step from the scalar result.
- [ ] **Lyapunov exponent relationship**: How do Lyapunov exponents of regular and $\alpha$-transformed maps relate? Run §4 of inverse_correspondence and characterize the $(\alpha, a)$ phase diagram.
- [ ] **Piecewise $\alpha(\omega)$**: Can we overcome the bootstrap failure of adaptive $\alpha$ by using different constant $\alpha$ values in different frequency windows?
- [ ] **Minimal discretization for RP resonances**: What resolution is needed for exact sub-leading Ruelle-Pollicott eigenvalues? Consider Jacobi polynomial basis weighted by arcsine measure (see [cross-theory connections §10](notebooks/cross_theory_connections.md)).

### Medium-term (requires new framework)

- [ ] **$\alpha(x)$ perturbation detection in experiment**: Apply the reconstruction formula $V'_\text{pert} = kx(\alpha(x) - 1)$ to real experimental data (e.g., optical trap perturbation measurement)
- [ ] **Chaos control comparison**: Quantitative comparison of $\alpha(x)$ chaos suppression vs. OGY and Pyragas methods — compute basins of attraction, robustness to noise, control effort
- [ ] **Non-equilibrium Green's functions**: Alpha-relaxation for the Kadanoff-Baym equations (two-time Dyson equation on the Keldysh contour)

### Speculative (may or may not lead anywhere)

- [ ] **Gauge-theoretic formulation**: Does treating $\alpha(x)$ as a gauge field (with the operator group $\Gamma_\alpha \Gamma_\beta = \Gamma_{\alpha\beta}$) produce conservation laws or topological invariants?
- [ ] **Concavity inversion in optimization**: Since $g'' = \alpha f''$, negative $\alpha$ flips convexity. Can this be used to escape local minima in non-convex optimization by temporarily inverting the landscape?
- [ ] **Koopman operator connection**: Does the $\alpha$-transform correspond to a similarity transformation on the Koopman operator? If so, the spectral invariance (equal Lyapunov exponents) follows immediately.
- [ ] **Renormalization group connection**: The Migdal-Kadanoff map is already an RG transformation. What does $\alpha$-transforming an RG map mean physically?

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
