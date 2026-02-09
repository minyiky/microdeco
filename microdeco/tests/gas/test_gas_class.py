"""
Tests for the Gas class in microdeco.gas module.
"""

import pytest

from microdeco.test_utils import approx
from microdeco.gas import Gas


# Test cases for Gas class
@pytest.mark.parametrize(
    "p_o2, p_n2, p_he, expected_name",
    [
        (21, 79, 0, "Air"),
        (100, 0, 0, "Oxygen"),
        (32, 68, 0, "EAN32"),
        (40, 60, 0, "EAN40"),
        (21, 0, 79, "Heliox21"),
        (21, 35, 44, "Tri21/44"),
    ],
)
def test_gas_class(
    p_o2: int,
    p_n2: float,
    p_he: float,
    expected_name: str,
):
    """
    Test the behaviour of the Gas class using precomputed values.
    """
    gas = Gas(p_o2=p_o2, p_n2=p_n2, p_he=p_he)
    assert gas.name == expected_name
    assert approx(gas.p_o2, p_o2 / 100.0)
    assert approx(gas.p_n2, p_n2 / 100.0)
    assert approx(gas.p_he, p_he / 100.0)


# Test the default gas is Air
def test_gas_default():
    """
    Test that the default gas is Air.
    """
    gas = Gas()
    assert gas.name == "Air"
    assert approx(gas.p_o2, 0.21)
    assert approx(gas.p_n2, 0.79)
    assert approx(gas.p_he, 0.0)


# Test invalid gas composition
def test_gas_invalid_composition():
    """
    Test that the Gas class raises an assertion error for invalid compositions.
    """
    with pytest.raises(AssertionError):
        Gas(p_o2=50, p_n2=30, p_he=30)  # Sum is 110%
    with pytest.raises(AssertionError):
        Gas(p_o2=50, p_n2=60, p_he=-10)  # Negative percentage
