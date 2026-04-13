"""
Alpha-transform and related stabilization methods.

The alpha-transform g(x) = alpha*f(x) + (1 - alpha)*x blends a map f
with the identity, converting unstable fixed points to stable ones
when alpha is chosen appropriately.
"""

import numpy as np
import sympy as sym


def alpha_transform(f, alpha, a, x):
    """
    Apply the alpha-transform: g(x) = alpha*f(a,x) + (1-alpha)*x.

    Parameters
    ----------
    alpha : float or callable
        If callable, alpha(x) is evaluated at each x (position-dependent
        relaxation). If float, standard constant-alpha transform.
    """
    if callable(alpha):
        a_val = alpha(x)
    else:
        a_val = alpha
    return a_val * f(a, x) + (1.0 - a_val) * x


def compute_optimal_alpha(f_symbolic, a_val, x_domain=(-2.0, 5.0), boundary=False):
    """
    Compute the optimal alpha for stabilizing unstable fixed points.

    Given f(a, x), finds the extremal slopes (via second derivative roots)
    and computes:
        alpha = 1 / (1 - df_max)   if df_max > 1   (optimal, zeroes g')
        alpha = 1 / (1 - df_min)   if df_min < -1   (optimal, zeroes g')

    With boundary=True, uses 2/(1-df) instead (marginal stability, g'=-1).

    Parameters
    ----------
    f_symbolic : callable
        A function f(a, x) that works with sympy symbols.
    a_val : float
        Numeric value of the parameter a.
    x_domain : tuple
        Domain (x_min, x_max) for numeric fallback when f'' has no roots.
    boundary : bool
        If True, return the boundary alpha (g'=-1) instead of optimal (g'=0).

    Returns
    -------
    alpha : float
        Optimal transform parameter.
    info : dict
        Dictionary with keys 'df_max', 'df_min', 'critical_points', 'method'.
    """
    a = sym.Symbol('a')
    x = sym.Symbol('x')

    f_expr = f_symbolic(a, x)
    df = sym.diff(f_expr, x)
    ddf = sym.diff(df, x)

    # Substitute numeric a
    df_num = sym.simplify(df.subs(a, a_val))
    ddf_num = sym.simplify(ddf.subs(a, a_val))

    # Find critical points of derivative (roots of second derivative)
    critical_x = sym.solve(ddf_num, x)

    df_max = 0.0
    df_min = 0.0
    slopes = []

    for xc in critical_x:
        slope = float(df_num.subs(x, xc))
        slopes.append((float(xc), slope))
        if np.isfinite(slope):
            df_max = max(df_max, slope)
            df_min = min(df_min, slope)

    # Fallback: if no critical points found (e.g. linear f'), evaluate
    # f' on a grid to find extremal slopes numerically
    if len(critical_x) == 0:
        x_grid = np.linspace(x_domain[0], x_domain[1], 2000)
        df_func = sym.lambdify(x, df_num, 'numpy')
        df_vals = df_func(x_grid)
        df_max = float(np.max(df_vals))
        df_min = float(np.min(df_vals))
        slopes = [('numeric_scan', df_max), ('numeric_scan', df_min)]

    coeff = 2.0 if boundary else 1.0

    info = {
        'df_max': df_max,
        'df_min': df_min,
        'critical_points': slopes,
        'method': None,
        'boundary': boundary,
    }

    if df_max > 1.0:
        info['method'] = 'max'
        return coeff / (1.0 - df_max), info
    elif df_min < -1.0:
        info['method'] = 'min'
        return coeff / (1.0 - df_min), info
    else:
        info['method'] = 'none_needed'
        return 1.0, info


def alpha_range_for_logistic(a):
    """
    For the logistic map f(x) = ax(1-x) with f'(x) = a - 2ax,
    the valid alpha range is (2/(1-a), 2/(1+a)) when a > 1.
    """
    if a <= 1:
        return None, None
    return 2.0 / (1.0 - a), 2.0 / (1.0 + a)


def optimal_alpha_at_x(f_prime, a, x):
    """
    Compute the optimal alpha at position x given local slope f'(a, x).

    Returns alpha(x) = 1 / (1 - f'(a, x)), which zeroes g'(x) at x.
    Falls back to a bounded value when |f'| ≈ 1.

    Parameters
    ----------
    f_prime : callable
        Derivative f'(a, x).
    a : float
        Parameter value.
    x : float or ndarray
        Position(s).

    Returns
    -------
    alpha : float or ndarray
    """
    slope = f_prime(a, x)
    denom = 1.0 - slope
    # Avoid division by zero near slope = 1 (neutral fixed points)
    return np.where(np.abs(denom) > 1e-12, 1.0 / denom, 0.0)


def make_alpha_func(f_prime, a, smooth=True, cap=5.0):
    """
    Construct a callable alpha(x) from the local slope f'(a, x).

    The idea: at each x, choose alpha(x) = 1/(1 - f'(a,x)) to zero
    the transformed derivative locally. This overcomes the obstruction
    of Proposition 4: different fixed points get different alpha values
    automatically.

    Parameters
    ----------
    f_prime : callable
        Derivative f'(a, x) of the map.
    a : float
        Map parameter.
    smooth : bool
        If True, clip alpha to [-cap, cap] to prevent divergence
        near neutral fixed points (f' ≈ 1).
    cap : float
        Maximum |alpha| when smooth=True.

    Returns
    -------
    alpha_func : callable
        A function alpha(x) ready to pass to alpha_transform.
    """
    def alpha_func(x):
        slope = f_prime(a, np.asarray(x, dtype=float))
        denom = 1.0 - slope
        alpha = np.where(np.abs(denom) > 1e-12, 1.0 / denom, 0.0)
        if smooth:
            alpha = np.clip(alpha, -cap, cap)
        return float(alpha) if np.ndim(alpha) == 0 else alpha

    return alpha_func
