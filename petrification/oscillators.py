"""
Damped harmonic oscillator systems for perturbation analysis.

Provides simulation of particle dynamics in arbitrary potentials,
fitting of damped sinusoidal trajectories, and construction of
relaxation iteration maps for alpha-transform perturbation detection.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
#  Potential / force builders
# ---------------------------------------------------------------------------

def hooke_force(k, x_eq=0.0):
    """
    Hooke force: F(x) = -k·(x - x_eq).

    Parameters
    ----------
    k : float
        Spring constant.
    x_eq : float
        Equilibrium position (default 0).

    Returns
    -------
    F : callable
        Force function F(x).
    """
    def F(x):
        return -k * (np.asarray(x, dtype=float) - x_eq)
    return F


def hooke_potential(k, x, x_eq=0.0):
    """Hooke potential: V(x) = ½k·(x - x_eq)²."""
    return 0.5 * k * (np.asarray(x, dtype=float) - x_eq)**2


def quartic_force(epsilon):
    """Force from quartic potential V = ε·x⁴/4: F(x) = -ε·x³."""
    def F(x):
        return -epsilon * np.asarray(x, dtype=float)**3
    return F


def quartic_potential(epsilon, x):
    """Quartic potential: V(x) = ε·x⁴/4."""
    return epsilon * np.asarray(x, dtype=float)**4 / 4


def combine_forces(*forces):
    """Sum multiple force functions into a single callable."""
    def F(x):
        return sum(f(x) for f in forces)
    return F


# ---------------------------------------------------------------------------
#  ODE simulation
# ---------------------------------------------------------------------------

def simulate(V_prime, gamma, m, x0, v0, t_span, dt=0.001):
    """
    Simulate damped motion in an arbitrary potential.

    Solves  m·ẍ + γ·ẋ + V'(x) = 0  with RK45.

    Parameters
    ----------
    V_prime : callable
        Potential gradient V'(x).  For Hooke's law V'(x) = k·x.
    gamma : float
        Damping coefficient.
    m : float
        Mass.
    x0, v0 : float
        Initial position and velocity.
    t_span : tuple
        (t_start, t_end).
    dt : float
        Output time step.

    Returns
    -------
    t : ndarray
        Time points.
    x : ndarray
        Position.
    v : ndarray
        Velocity.
    """
    def ode(_t, y):
        pos, vel = y
        return [vel, -V_prime(pos) / m - (gamma / m) * vel]

    t_eval = np.arange(t_span[0], t_span[1], dt)
    sol = solve_ivp(ode, t_span, [x0, v0], t_eval=t_eval,
                    method='RK45', rtol=1e-12, atol=1e-14)
    return sol.t, sol.y[0], sol.y[1]


# ---------------------------------------------------------------------------
#  Trajectory analysis
# ---------------------------------------------------------------------------

def find_peaks(x):
    """Return indices of local maxima in a 1-D signal."""
    idx = []
    for i in range(1, len(x) - 1):
        if x[i] > x[i - 1] and x[i] > x[i + 1]:
            idx.append(i)
    return np.array(idx)


def fit_damped_sinusoid(t, x):
    """
    Fit  x(t) = A·exp(-β·t)·cos(ωd·t + φ) + x_eq  to data.

    Parameters
    ----------
    t : ndarray
        Time points.
    x : ndarray
        Position data.

    Returns
    -------
    params : dict
        Keys: A, beta, omega_d, phi, x_eq.
    x_fit : ndarray
        Fitted curve evaluated on t.
    """
    x_eq_g = np.mean(x[-max(1, len(x) // 10):])
    xc = x - x_eq_g
    A_g = float(np.max(np.abs(xc)))

    crossings = np.where(np.diff(np.sign(xc)))[0]
    if len(crossings) >= 2:
        omega_g = float(np.pi / np.median(np.diff(t[crossings])))
    else:
        omega_g = 1.0

    pk = find_peaks(x)
    if len(pk) >= 2:
        pv = np.abs(xc[pk])
        pos = pv > 1e-12
        if np.sum(pos) >= 2:
            pk_p, pv_p = pk[pos], pv[pos]
            beta_g = float(np.log(pv_p[0] / pv_p[-1])
                           / (t[pk_p[-1]] - t[pk_p[0]]))
            beta_g = max(beta_g, 0.01)
        else:
            beta_g = 0.1
    else:
        beta_g = 0.1

    def model(t, A, beta, omega_d, phi, x_eq):
        return A * np.exp(-beta * t) * np.cos(omega_d * t + phi) + x_eq

    try:
        popt, _ = curve_fit(
            model, t, x,
            p0=[A_g, beta_g, omega_g, 0.0, x_eq_g],
            maxfev=20000,
        )
    except RuntimeError:
        popt = [A_g, beta_g, omega_g, 0.0, x_eq_g]

    params = dict(zip(['A', 'beta', 'omega_d', 'phi', 'x_eq'], popt))
    x_fit = model(t, *popt)
    return params, x_fit


def extract_spring_constant(omega_d, beta, m):
    """
    Recover spring constant from fitted parameters.

    ωd² = k/m - β²  ⟹  k = m·(ωd² + β²).
    """
    return m * (omega_d**2 + beta**2)


# ---------------------------------------------------------------------------
#  Relaxation maps for fixed-point analysis
# ---------------------------------------------------------------------------

def relaxation_map(V_prime, eta):
    """
    Construct the relaxation iteration map  f(x) = x - η·V'(x).

    Fixed points satisfy V'(x) = 0.  The slope at a fixed point x*
    is f'(x*) = 1 - η·V''(x*).

    Parameters
    ----------
    V_prime : callable
        Potential gradient.
    eta : float
        Step size.

    Returns
    -------
    f : callable
        Iteration map, accepts scalar or ndarray.
    """
    def f(x):
        return x - eta * V_prime(np.asarray(x, dtype=float))
    return f


def measure_alpha_profile(f_model, f_true, x_range, n_points=500):
    """
    Measure position-dependent α(x) that transforms f_model into f_true.

    From  g(x) = α(x)·f_model(x) + (1 - α(x))·x = f_true(x),  solve::

        α(x) = [f_true(x) - x] / [f_model(x) - x]

    Parameters
    ----------
    f_model, f_true : callable
        Model and true iteration maps.
    x_range : tuple
        (x_min, x_max).
    n_points : int
        Number of evaluation points.

    Returns
    -------
    x_eval : ndarray
        Evaluation positions.
    alpha : ndarray
        α(x) at each position (NaN where f_model(x) ≈ x).
    """
    x_eval = np.linspace(x_range[0], x_range[1], n_points)
    h_model = f_model(x_eval) - x_eval
    h_true = f_true(x_eval) - x_eval
    with np.errstate(divide='ignore', invalid='ignore'):
        alpha = np.where(np.abs(h_model) > 1e-15, h_true / h_model, np.nan)
    return x_eval, alpha


def infer_perturbation(k_model, x_eval, alpha_profile):
    """
    Reconstruct the perturbation potential from the α(x) profile.

    From  α(x) = 1 + V'_pert(x)/(k_model·x)  we get::

        V'_pert(x) = k_model · x · (α(x) - 1)
        V_pert(x)  = ∫ V'_pert dx   (cumulative trapezoid)

    Parameters
    ----------
    k_model : float
        Known model spring constant.
    x_eval : ndarray
        Positions where α was measured.
    alpha_profile : ndarray
        Measured α(x).

    Returns
    -------
    V_prime_pert : ndarray
        Inferred perturbation gradient.
    V_pert : ndarray
        Inferred perturbation potential (zeroed near x = 0).
    """
    V_prime_pert = k_model * x_eval * (alpha_profile - 1)
    dx = np.gradient(x_eval)
    V_pert = np.cumsum(V_prime_pert * dx)
    i_zero = np.argmin(np.abs(x_eval))
    V_pert -= V_pert[i_zero]
    return V_prime_pert, V_pert
