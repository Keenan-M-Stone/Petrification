"""
Fixed-point iteration engines for discrete dynamical systems.

Provides the core iteration loop used for cobweb diagrams,
convergence analysis, and alpha-transform visualization.
"""

import numpy as np


def iterate(f, a, x0, n_iter=100):
    """
    Iterate x_{n+1} = f(a, x_n) and return the trajectory.

    Parameters
    ----------
    f : callable
        Map function f(a, x).
    a : float
        Parameter value.
    x0 : float
        Initial condition.
    n_iter : int
        Number of iterations.

    Returns
    -------
    trajectory : ndarray of shape (n_iter + 1,)
    """
    traj = np.zeros(n_iter + 1)
    traj[0] = x0
    x = x0
    for i in range(n_iter):
        x = f(a, x)
        traj[i + 1] = x
    return traj


def iterate_transformed(f, alpha, a, x0, n_iter=100):
    """
    Iterate the alpha-transformed map g(x) = alpha*f(a,x) + (1-alpha)*x.

    Parameters
    ----------
    alpha : float or callable
        If callable, alpha(x) is evaluated at each step.
    """
    traj = np.zeros(n_iter + 1)
    traj[0] = x0
    x = x0
    for i in range(n_iter):
        a_val = alpha(x) if callable(alpha) else alpha
        x = a_val * f(a, x) + (1.0 - a_val) * x
        traj[i + 1] = x
    return traj


def cobweb_data(f, a, x0, n_iter=20, alpha=None):
    """
    Generate (x, y) data for a cobweb diagram.

    Returns arrays Ix, Iy that when plotted show the staircase
    pattern of fixed-point iteration: start at (x0, 0), go vertical
    to (x0, f(x0)), horizontal to (f(x0), f(x0)), etc.

    Parameters
    ----------
    f : callable
        Map function f(a, x).
    a : float
        Parameter value.
    x0 : float
        Initial condition.
    n_iter : int
        Number of iteration steps.
    alpha : float or None
        If provided, use the alpha-transform.

    Returns
    -------
    Ix, Iy : ndarray
        x and y coordinates for the cobweb plot.
    """
    def step(x):
        if alpha is not None:
            a_val = alpha(x) if callable(alpha) else alpha
            return a_val * f(a, x) + (1.0 - a_val) * x
        return f(a, x)

    Ix = [x0, x0]
    Iy = [0.0, step(x0)]

    x = step(x0)
    for _ in range(n_iter):
        x_new = step(x)
        Ix.extend([x, x])
        Iy.extend([x, x_new])
        x = x_new

    return np.array(Ix), np.array(Iy)


def find_fixed_points(f, a, x_range, n_points=1000, tol=1e-10):
    """
    Numerically find fixed points of f(a, x) = x in a given range.

    Returns approximate fixed-point locations.
    """
    x = np.linspace(x_range[0], x_range[1], n_points)
    g = f(a, x) - x
    # Find sign changes
    sign_changes = np.where(np.diff(np.sign(g)))[0]

    fps = []
    for idx in sign_changes:
        # Linear interpolation for better estimate
        x1, x2 = x[idx], x[idx + 1]
        g1, g2 = g[idx], g[idx + 1]
        if abs(g2 - g1) > 0:
            xfp = x1 - g1 * (x2 - x1) / (g2 - g1)
            fps.append(xfp)

    return np.array(fps)
