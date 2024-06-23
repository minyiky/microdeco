"""
This module provides functions for calculating tissue
pressure using various equations.
"""

from .gas_equations import schreiner, buhlmann

__all__ = ["schreiner", "buhlmann"]
