# Citations — Petrification Project

**Compiled April 2026**

All references cited in Petrification notebooks, scripts, and documentation.
Only references explicitly appearing in the source files are included.

---

## 1. Core Theory & Fixed-Point Iteration

### Krasnoselskii-Mann Relaxation (foundation of alpha-transform equivalence)

- **Krasnoselskii, M.A.** "Two remarks on the method of successive approximations." *Uspekhi Mat. Nauk* 10(1), 123–127 (1955).
  - *Cited in:* alpha_transform.ipynb
  - The alpha-transform $g(x) = \alpha f(x) + (1-\alpha)x$ is the Krasnoselskii-Mann iteration family.

### Textbooks

- **Strogatz, S.H.** *Nonlinear Dynamics and Chaos.* 2nd ed., CRC Press, 2015.
  - *Cited in:* alpha_transform.ipynb, eigenstate_fixedpoint_correspondence.ipynb
  - Bifurcation theory, period-doubling, cobweb diagrams.

- **Trefethen, L.N. and Bau, D.** *Numerical Linear Algebra.* SIAM, 1997.
  - *Cited in:* alpha_transform.ipynb, eigenstate_fixedpoint_correspondence.ipynb, turbiner_riccati_bloch.ipynb
  - Power iteration, eigenvalue computation, matrix methods.

### Nonlinear Eigenvalue Problems

- **Güttel, S. and Tisseur, F.** "The nonlinear eigenvalue problem." *Acta Numerica* 26, 1–94 (2017).
  - *Cited in:* alpha_transform.ipynb

---

## 2. Quantum Mechanics & Spectral Theory

### Turbiner's Nonlinearization (Riccati-Bloch)

- **Turbiner, A.V.** "The eigenvalue spectrum in quantum mechanics and the nonlinearization procedure." *Soviet Physics Uspekhi* 27(9), 668–694 (1984).
  - *Cited in:* eigenstate_fixedpoint_correspondence.ipynb
  - Riccati-Bloch equation, eigenvalue detection via bounded trajectories.

- **Turbiner, A.V.** "Quasi-exactly-solvable problems and sl(2) algebra." *Communications in Mathematical Physics* 118(3), 467–474 (1988).
  - *Cited in:* turbiner_riccati_bloch.ipynb, eigenstate_fixedpoint_correspondence.ipynb

- **Turbiner, A.V.** "One-dimensional quasi-exactly solvable Schrödinger equations." *Physics Reports* 642, 1–71 (2016).
  - *Cited in:* turbiner_riccati_bloch.ipynb

- **Turbiner, A.V. and Ushveridze, A.G.** "Anharmonic oscillator: constructing the strong coupling expansions." *Journal of Mathematical Physics* 29, 2053–2062 (1988).
  - *Cited in:* turbiner_riccati_bloch.ipynb

### Quantum Mechanics Textbooks

- **Griffiths, D.J. and Schroeter, D.F.** *Introduction to Quantum Mechanics.* 3rd ed., Cambridge University Press, 2018.
  - *Cited in:* multiple notebooks (general reference)

- **Bender, C.M. and Orszag, S.A.** *Advanced Mathematical Methods for Scientists and Engineers.* Springer, 1999.
  - *Cited in:* turbiner_riccati_bloch.ipynb
  - WKB approximation, semiclassical quantization, PT-symmetric quantum mechanics.

### Riccati Equations

- **Reid, W.T.** *Riccati Differential Equations.* Academic Press, 1972.
  - *Cited in:* alpha_transform.ipynb, eigenstate_fixedpoint_correspondence.ipynb, turbiner_riccati_bloch.ipynb

---

## 3. Dynamical Systems & Chaos Control

### OGY Method

- **Ott, E., Grebogi, C., and Yorke, J.A.** "Controlling chaos." *Physical Review Letters* 64, 1196–1199 (1990).
  - *Cited in:* inverse_correspondence.ipynb, chaos_control_comparison.ipynb

### Pyragas Delayed Feedback Control

- **Pyragas, K.** "Continuous control of chaos by self-controlling feedback." *Physics Letters A* 170, 421–428 (1992).
  - *Cited in:* inverse_correspondence.ipynb, chaos_control_comparison.ipynb

### Odd Number Limitation

- **Nakajima, H.** "On analytical properties of delayed feedback control of chaos." *Physics Letters A* 232, 207–210 (1997).
  - *Cited in:* chaos_control_comparison.ipynb
  - Pyragas DFC with positive K cannot stabilize orbits with odd number of real Floquet multipliers > 1.

- **Just, W., Bernard, T., Ostheimer, M., Reibold, E., and Benner, H.** "Mechanism of time-delayed feedback control." *Physical Review Letters* 78, 203–206 (1997).
  - *Cited in:* chaos_control_comparison.ipynb

---

## 4. Transfer Operators & Spectral Analysis

### Koopman Operator Theory

- **Koopman, B.O.** "Hamiltonian systems and transformation in Hilbert space." *Proceedings of the National Academy of Sciences* 17(5), 315–318 (1931).
  - *Cited in:* eigenstate_fixedpoint_correspondence.ipynb

- **Mezić, I.** "Spectral properties of dynamical systems, model reduction and decompositions." *Nonlinear Dynamics* 41, 309–325 (2005).
  - *Cited in:* eigenstate_fixedpoint_correspondence.ipynb, turbiner_riccati_bloch.ipynb

- **Budišić, M., Mohr, R., and Mezić, I.** "Applied Koopmanism." *Chaos* 22(4), 047510 (2012).
  - *Cited in:* eigenstate_fixedpoint_correspondence.ipynb

### Ruelle-Pollicott Resonances

- **Ruelle, D.** "Resonances of chaotic dynamical systems." *Physical Review Letters* 56, 405–407 (1986).
  - *Cited in:* turbiner_riccati_bloch.ipynb

- **Ruelle, D.** "Zeta-functions for expanding maps and Anosov flows." *Inventiones Mathematicae* 34, 231–242 (1976).
  - *Cited in:* eigenstate_fixedpoint_correspondence.ipynb

### Transfer Operator Approximation

- **Froyland, G.** "Approximating physical invariant measures of mixing dynamical systems." *Nonlinearity* 11, 1043–1090 (1998).
  - *Cited in:* turbiner_riccati_bloch.ipynb

- **Blank, M., Keller, G., and Liverani, C.** "Ruelle-Perron-Frobenius spectrum for Anosov maps." *Nonlinearity* 15, 1905–1919 (2002).
  - *Cited in:* turbiner_riccati_bloch.ipynb, eigenstate_fixedpoint_correspondence.ipynb

---

## 5. Statistical Mechanics & Renormalization Group

- **Wilson, K.G.** "The renormalization group: Critical phenomena and the Kondo problem." *Reviews of Modern Physics* 47(4), 773–840 (1975).
  - *Cited in:* alpha_transform.ipynb, eigenstate_fixedpoint_correspondence.ipynb

- **Kadanoff, L.P.** "Scaling laws for Ising models near $T_c$." *Physics Physique Fizika* 2(6), 263–283 (1966).
  - *Cited in:* alpha_transform.ipynb, eigenstate_fixedpoint_correspondence.ipynb

- **Migdal, A.A.** "Recursion equations in gauge theories." *Soviet Physics JETP* 42(3), 413–418 (1975).
  - *Cited in:* alpha_transform.ipynb, eigenstate_fixedpoint_correspondence.ipynb

---

## 6. Self-Consistent Field Methods (Dyson Equation Context)

These references establish that empirical mixing parameters are standard practice but lack theoretical derivation — the gap that the crossover notebook's $\alpha^* \approx 0.5$ result addresses.

- **Chibani, S. et al.** "Self-consistent Green's function embedding for advanced electronic structure methods." *Physical Review B* (2016).
  - *Cited in:* alpha_transform.ipynb, crossover_alpha_turbiner.ipynb (context)
  - Uses empirical mixing $\lambda = 0.5$ in GW calculations.

- **Caruso, F. et al.** *Physical Review B* (2013).
  - *Cited in:* alpha_transform.ipynb
  - Adaptive mixing in GW calculations.

- **Forster.** *Journal of Chemical Theory and Computation* (2025).
  - *Cited in:* alpha_transform.ipynb
  - Adaptive mixing methods.

- **Tadano, T. and Tsuneyuki, S.** *Physical Review B* (2015).
  - *Cited in:* alpha_transform.ipynb
  - Uses $\alpha = 0.1$ for phonon self-consistency.

---

## 7. Data-Driven Methods

- **Brunton, S.L. and Kutz, J.N.** *Data-Driven Science and Engineering.* 2nd ed., Cambridge University Press, 2022.
  - *Cited in:* eigenstate_fixedpoint_correspondence.ipynb

---

## 8. Chaos Control

Used in `chaos_control_comparison.ipynb` (OGY method, Pyragas feedback, odd-number limitation):

- **Ott, E., Grebogi, C., and Yorke, J.A.** "Controlling chaos." *Physical Review Letters* **64**, 1196–1199 (1990). DOI:10.1103/PhysRevLett.64.1196  
  *Foundational OGY method — small parameter perturbations to stabilise embedded unstable periodic orbits.*

- **Pyragas, K.** "Continuous control of chaos by self-controlling feedback." *Physics Letters A* **170**, 421–428 (1992). DOI:10.1016/0375-9601(92)90745-8  
  *Delayed feedback control of chaos (TDAS); the Pyragas method benchmarked against α(x) control.*

- **Nakajima, H.** "On analytical properties of delayed feedback control of chaos." *Physics Letters A* **232**, 207–210 (1997). DOI:10.1016/S0375-9601(97)00362-9  
  *Proved the odd-number limitation: delayed feedback cannot stabilise orbits with an odd number of real Floquet multipliers greater than 1. Explains the 0% basin in our Pyragas experiments.*

- **Just, W., Bernard, T., Ostheimer, M., Reibold, E., and Benner, H.** "Mechanism of time-delayed feedback control." *Physical Review Letters* **78**, 203–206 (1997). DOI:10.1103/PhysRevLett.78.203  
  *Provides the mechanism underlying the odd-number limitation; corroborates Nakajima (1997).*

---

## 9. Ergodic Theory / Entropy

Used implicitly across `lyapunov_alpha_relationship.ipynb`, `rp_resonance_basis.ipynb`:

- **Pesin, Ya.B.** "Characteristic Lyapunov exponents and smooth ergodic theory." *Russian Mathematical Surveys* **32**, 55–114 (1977). DOI:10.1070/RM1977v032n04ABEH001639  
  *Pesin's theorem: for smooth ergodic systems, KS (metric) entropy equals the sum of positive Lyapunov exponents. Basis for the Λ = KS entropy statement in our Lyapunov notebook.*

---

## 10. Related Work on α-Parametrised Maps (for Distinction)

- **Pires, M.A., Tsallis, C., and Curado, E.M.F.** "Composing α-Gauss and logistic maps: Gradual and sudden transitions to chaos." *Physical Review E* **112**, 034209 (2025). DOI:10.1103/lwfn-qrjt  
  *Introduces the α-Gauss-Logistic map $x_{t+1} = f_L(x_t)x_t^{-\alpha} - \lfloor f_L(x_t)x_t^{-\alpha}\rfloor$. Different from our α-transform (which is the Krasnoselskii-Mann convex combination $g_\alpha = \alpha f + (1-\alpha)\text{id}$): their α parametrises a Gauss-map composition, not a relaxation step. Cited for context; no conflict with our work.*

---

## 11. Cross-Theory Connections (from cross_theory_connections.md)

Frameworks noted as connections but with incomplete bibliographic references:

- **Kottos, T. and Smilansky, U.** Quantum chaos on graphs (1999).
- **Boyland / Hall** — Braid theory of periodic orbits (1990s).
- **Strogatz, S.H.** *Nonlinear Dynamics and Chaos*, 2nd ed. (2015). §10.4 gives the closed-form logistic orbit for $a=4$ via the $x=\sin^2(\pi\theta)$ conjugacy.

---

## Summary

| Category | Count |
|----------|-------|
| Core theory / iteration | 4 |
| Quantum mechanics | 7 |
| Chaos control | 4 (expanded to 8 with Nakajima/Just) |
| Transfer operators | 6 |
| Statistical mechanics / RG | 3 |
| Self-consistent field methods | 4 |
| Data-driven methods | 1 |
| Ergodic theory | 1 |
| Related α-map work (for distinction) | 1 |
| Incomplete references | 3 |
| **Total** | **≥ 39** |

Most cited author: **Turbiner, A.V.** (4 distinct works).
Time span: 1931–2025.
Primary disciplines: dynamical systems, quantum mechanics, numerical analysis, chaos theory.
