"""
Bifurcation diagram generation.

Computes long-term iterates of a map f(a,x) over a range of parameter
values to produce bifurcation diagrams, with optional alpha-transform.
"""

import numpy as np


def compute_bifurcation(f, a_range, x0_samples, n_iter=500,
                        n_discard=100, x_max=5.0):
    """
    Compute bifurcation data for map f(a, x).

    Parameters
    ----------
    f : callable
        Map function f(a, x).
    a_range : array-like
        Parameter values to sweep.
    x0_samples : array-like
        Initial conditions to try for each a.
    n_iter : int
        Total iterations per (a, x0) pair.
    n_discard : int
        Transient iterations to discard before recording.
    x_max : float
        Divergence threshold.

    Returns
    -------
    a_vals : ndarray
        Parameter values (repeated for each attractor point).
    x_vals : ndarray
        Corresponding long-term iterates.
    """
    a_out = []
    x_out = []

    for a in a_range:
        attractor = set()
        for x0 in x0_samples:
            x = x0
            diverged = False
            for i in range(n_iter):
                x = f(a, x)
                if abs(x) > x_max:
                    diverged = True
                    break
                if i >= n_discard:
                    # Round to avoid float noise in the set
                    attractor.add(round(x, 10))

            if diverged:
                continue

        for xv in attractor:
            a_out.append(a)
            x_out.append(xv)

    return np.array(a_out), np.array(x_out)


def compute_bifurcation_transformed(f, alpha, a_range, x0_samples,
                                    n_iter=500, n_discard=100, x_max=5.0):
    """
    Compute bifurcation data for the alpha-transformed map
    g(x) = alpha*f(a,x) + (1-alpha)*x.

    Parameters are the same as compute_bifurcation, plus alpha.
    """
    def g(a, x):
        return alpha * f(a, x) + (1.0 - alpha) * x

    return compute_bifurcation(g, a_range, x0_samples,
                               n_iter, n_discard, x_max)
