"""
Equations used to control gas loading simultaion in tissues.

This file covers
    - The Schreiner equation
    - The Buhlmann equation
"""

import math


def schreiner(
    p_alv: float,
    p_i: float,
    k: float,
    t: float,
    R: float,  # pylint: disable=invalid-name
) -> float:
    """
    Calculate changes in tissue pressure using the Schreiner equation.

    This function calculates the Schreiner equation using the provided
    parameters.

    Parameters:
        p_alv (float): The absolute partial pressure of inspired gas.
        p_i (float): The initial absolute partial pressure of the
            gas in the tissue.
        k (float): Gas decay constant (tissue specific).
        t (float): The time in minutes.
        R (float): The rate of change of ambient pressure (bar per minute).

    Returns:
        float: The final pressure of the tissue.
    """
    exponent = math.exp(-(k * t))
    return p_alv + R * (t - 1 / k) - (p_alv - p_i - (R / k)) * exponent


def buhlmann(
    gf: float,
    p_n2: float,
    p_he: float,
    A_n2: float,  # pylint: disable=invalid-name
    A_he: float,  # pylint: disable=invalid-name
    B_n2: float,  # pylint: disable=invalid-name
    B_he: float,  # pylint: disable=invalid-name
) -> float:
    """
    Calculates pressure ceiling of a tissue using the Buhlmann equation.

    This uses gradient factors to proivde an additional safety margin.

    Args:
        gf (float): The gradient factor.
        p_n2 (float): The partial pressure of Nitrogen.
        p_he (float): The partial pressure of Helium.
        A_n2 (float): The A coefficient for Nitrogen.
        A_he (float): The A coefficient for Helium.
        B_n2 (float): The B coefficient for Nitrogen.
        B_he (float): The B coefficient for Helium.

    Returns:
        float: The safe pressure ceiling of the tissue.

    Raises:
        None
    """
    p = p_n2 + p_he
    A = (A_n2 * p_n2 + A_he * p_he) / p  # pylint: disable=invalid-name
    B = (B_n2 * p_n2 + B_he * p_he) / p  # pylint: disable=invalid-name
    return (p - A * gf) / ((gf / B) + 1 - gf)
