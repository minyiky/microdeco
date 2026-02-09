"""
Tests for the GradientFactor class in microdeco.models module.
"""

import pytest

from microdeco.models import GradientFactors


# Error cases for GradientFactor class
@pytest.mark.parametrize(
    "gf_low, gf_high",
    [
        (-0.5, 0.5),  # gf_low too low
        (0.5, 1.5),  # gf_high too high
        (0.7, 0.5),  # gf_low greater than gf_high
    ],
    ids=[
        "gf_low < 0",
        "gf_high > 1",
        "gf_low > gf_high",
    ],
)
def test_gradient_factor_invalid(gf_low: float, gf_high: float):
    """
    Test that the GradientFactor class raises a ValueError for invalid inputs.
    """
    with pytest.raises(ValueError):
        GradientFactors(gf_low=gf_low, gf_high=gf_high)
