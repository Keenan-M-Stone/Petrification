"""
Quantum potential functions for eigenvalue problems.
"""

import numpy as np


def harmonic(x, omega=1.0):
    """Harmonic oscillator: V(x) = (1/2)*omega^2*x^2."""
    return 0.5 * omega**2 * x**2


def anharmonic(x, g=0.1):
    """Quartic anharmonic oscillator: V(x) = x^2 + g*x^4."""
    return x**2 + g * x**4


def double_well(x, a=1.0):
    """Symmetric double well: V(x) = (x^2 - a^2)^2."""
    return (x**2 - a**2)**2


def finite_square_well(x, V0=50.0, a=1.0):
    """Finite square well: V=0 for |x|<a, V=V0 otherwise."""
    return np.where(np.abs(x) < a, 0.0, V0)


def morse(x, D=10.0, alpha=1.0, x0=0.0):
    """Morse potential: V(x) = D*(1 - exp(-alpha*(x-x0)))^2."""
    return D * (1.0 - np.exp(-alpha * (x - x0)))**2


def coulomb(r, Z=1.0, ell=0):
    """
    Effective radial Coulomb potential for the hydrogen-like atom.

    After substituting u(r) = r*psi(r), the radial Schrödinger equation
    becomes  -u'' + V_eff(r) u = E u  (in atomic units hbar=m=e=1) with:

        V_eff(r) = -Z/r + ell*(ell+1)/(2*r^2)

    Exact eigenvalues (for ell=0): E_n = -Z^2 / (2*n^2), n=1,2,...

    Parameters
    ----------
    r : ndarray
        Radial coordinate (must be > 0).
    Z : float
        Nuclear charge.
    ell : int
        Angular momentum quantum number.
    """
    r = np.asarray(r, dtype=float)
    V = np.where(r > 0, -Z / r + ell * (ell + 1) / (2.0 * r**2), 0.0)
    return V


# Registry of named potentials
POTENTIAL_REGISTRY = {
    "harmonic": harmonic,
    "anharmonic": anharmonic,
    "double_well": double_well,
    "finite_square_well": finite_square_well,
    "morse": morse,
    "coulomb": coulomb,
}


def get_potential(name):
    """Retrieve a potential function by name."""
    if name not in POTENTIAL_REGISTRY:
        available = ", ".join(POTENTIAL_REGISTRY.keys())
        raise KeyError(f"Unknown potential '{name}'. Available: {available}")
    return POTENTIAL_REGISTRY[name]
