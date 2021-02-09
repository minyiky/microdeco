from .engine import Engine

__version__ = '0.0.1'


def create(*args, **kwargs):
    """
    Create decompression engine .

    The decompression model validation is enabled by default.

    Usage

    >>> import decotengu
    >>> engine = decotengu.create()
    >>> engine.add_gas(0, 21)
    >>> data = list(engine.calculate(35, 40))
    >>> engine.deco_table.total
    44.0

    :param time_delta: Time between dive steps.
    :param validate: Validate decompression data with decompression model
                     validator.
    """
    engine = Engine(**kwargs)

    return engine


