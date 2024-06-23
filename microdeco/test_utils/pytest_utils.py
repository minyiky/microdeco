"""
Utitilies to help with using the pytest framework
"""

import pytest


def approx(result: float, *args, **kwargs) -> bool:  # type: ignore
    """
    Compare the given result with the expected value
    using pytest's `approx` function.

    Use to aviod tpye comparison errors.

    Args:
        result (float): The value to be compared.
        *args: Variable length argument list to be passed to `pytest.approx`.
        **kwargs: Keyword argument dictionary to be passed to `pytest.approx`.

    Returns:
        bool: True if the result is approximately equal to the expected value,
            False otherwise.
    """
    return result == pytest.approx(*args, **kwargs)  # type: ignore
