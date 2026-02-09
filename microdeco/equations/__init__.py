"""
This module provides functions for calculating tissue
pressure using various equations.
"""

from .gas_equations import (
    load_tissue,
    load_tissue_constant,
    buhlmann,
)

__all__ = ["load_tissue_constant", "load_tissue", "buhlmann"]
