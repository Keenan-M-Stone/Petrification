"""
Alpha-transform and related stabilization methods.

The alpha-transform g(x) = alpha*f(x) + (1 - alpha)*x blends a map f
with the identity, converting unstable fixed points to stable ones
when alpha is chosen appropriately.
"""

import numpy as np
import sympy as sym


def alpha_transform(f, alpha, a, x):
    """Apply the alpha-transform: g(x) = alpha*f(a,x) + (1-alpha)*x."""
    return alpha * f(a, x) + (1.0 - alpha) * x


def compute_optimal_alpha(f_symbolic, a_val, x_domain=(-2.0, 5.0)):
    """
    Compute the optimal alpha for stabilizing unstable fixed points.

    Given f(a, x), finds the extremal slopes (via second derivative roots)
    and computes:
        alpha = 2 / (1 - df_max)   if df_max > 1
        alpha = 2 / (1 - df_min)   if df_min < -1

    Parameters
    ----------
    f_symbolic : callable
        A function f(a, x) that works with sympy symbols.
    a_val : float
        Numeric value of the parameter a.
    x_domain : tuple
        Domain (x_min, x_max) for numeric fallback when f'' has no roots.

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

    info = {
        'df_max': df_max,
        'df_min': df_min,
        'critical_points': slopes,
        'method': None,
    }

    if df_max > 1.0:
        info['method'] = 'max'
        return 2.0 / (1.0 - df_max), info
    elif df_min < -1.0:
        info['method'] = 'min'
        return 2.0 / (1.0 - df_min), info
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
