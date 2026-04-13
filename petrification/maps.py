"""
Discrete maps for fixed-point analysis.

Contains the logistic map, Migdal-Kadanoff renormalization map,
and a registry for user-defined maps.
"""

import numpy as np


def logistic(a, x):
    """Standard logistic map: f(x) = a*x*(1-x)."""
    return a * x * (1.0 - x)


def migdal_kadanoff(a, x):
    """Migdal-Kadanoff renormalization map: f(x) = a*x^2 / (1+x^2)^2."""
    return a * x**2 / (1.0 + x**2)**2


def exponential_saturation(a, x):
    """Exponential saturation map: f(x) = a*(1 - exp(-x))."""
    return a * (1.0 - np.exp(-x))


# Registry of named maps for easy switching
MAP_REGISTRY = {
    "logistic": logistic,
    "migdal_kadanoff": migdal_kadanoff,
    "exponential_saturation": exponential_saturation,
}


def get_map(name):
    """Retrieve a map function by name from the registry."""
    if name not in MAP_REGISTRY:
        available = ", ".join(MAP_REGISTRY.keys())
        raise KeyError(f"Unknown map '{name}'. Available: {available}")
    return MAP_REGISTRY[name]


def register_map(name, func):
    """Register a custom map function. Must accept (a, x) arguments."""
    MAP_REGISTRY[name] = func
