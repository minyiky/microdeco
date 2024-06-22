import pytest

def approx(result: float, *args, **kwargs) -> bool: # type: ignore
    return result == pytest.approx(*args, **kwargs) # type: ignore