"""
Test cases for microdeco.equations
"""

import pytest

from microdeco.test_utils import approx
from microdeco.equations import schreiner, buhlmann


# Test cases for schreiner function
@pytest.mark.parametrize(
    "p_alv, p_t, k, t, R, expected",
    [
        (0.637364, 0.74065446, 0.138629, 1.5, 1.36, 0.919397),  # descent
        (2.677364, 0.919397, 0.138629, 20, 0, 2.567490),  # constant_depth
        (2.677364, 2.567490, 0.138629, 2, -0.68, 2.421830),  # ascent
    ],
)
def test_schreiner_function(
    p_alv: float,
    p_t: float,
    k: float,
    t: float,
    R: float,  # pylint: disable=invalid-name
    expected: float,
):
    """
    Test the behaviour of the schreiner function using precomputed values.
    """
    assert approx(
        schreiner(p_alv=p_alv, p_i=p_t, k=k, t=t, R=R),
        expected,
        rel=1e-5,
    )


# Test cases for buhlmann function
@pytest.mark.parametrize(
    "gf, p_n2, p_he, A_n2, A_he, B_n2, B_he, expected",
    [
        (0.3, 0.74065446, 0, 1.1696, 0, 0.5578, 0, 0.314886),  # surface
        (0.3, 0.919397, 0, 1.1696, 0, 0.5578, 0, 0.4592862),  # descent
        (0.3, 2.567490, 0, 1.1696, 0, 0.5578, 0, 1.7907266),  # constant depth
        (0.3, 2.421840, 0, 1.1696, 0, 0.5578, 0, 1.6730607),  # ascent
    ],
)
def test_buhlmann_function(
    gf: float,
    p_n2: float,
    p_he: float,
    A_n2: float,  # pylint: disable=invalid-name
    A_he: float,  # pylint: disable=invalid-name
    B_n2: float,  # pylint: disable=invalid-name
    B_he: float,  # pylint: disable=invalid-name
    expected: float,
):
    """
    Test the behaviour of the buhlmann function using precomputed values.
    """
    assert approx(
        buhlmann(gf, p_n2, p_he, A_n2, A_he, B_n2, B_he),
        expected,
        rel=1e-5,
    )
