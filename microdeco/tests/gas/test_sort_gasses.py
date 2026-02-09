"""
Tests for the sort_gasses function in microdeco.gas module.
"""

import pytest

from microdeco.test_utils import approx
from microdeco.gas import Gas
import itertools

AIR = Gas(21, 79, 0)
OXYGEN = Gas(100, 0, 0)
EAN32 = Gas(32, 68, 0)
EAN40 = Gas(40, 60, 0)
HELIOX21 = Gas(21, 0, 79)
HELIOX32 = Gas(32, 0, 68)
TRI21_35 = Gas(21, 44, 35)


# Test cases for sort_gasses function
@pytest.mark.parametrize(
    "expected_order",
    [
        ([OXYGEN, EAN40, EAN32, AIR]),
        ([OXYGEN, HELIOX32, HELIOX21]),
        ([HELIOX21, TRI21_35, AIR]),
        ([EAN32, TRI21_35, AIR]),
    ],
)
def test_sort_gasses(expected_order: list[Gas]):
    """
    Test the behaviour of the sort_gasses function using precomputed values.
    """
    for perm in itertools.permutations(expected_order):
        sorted_gases = sorted(list(perm))
        assert sorted_gases == expected_order
