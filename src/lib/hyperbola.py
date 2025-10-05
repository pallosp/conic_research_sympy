from collections.abc import Sequence

from sympy import Expr, Matrix

from lib.central_conic import ConicFromFociAndRadius
from lib.distance import PointPointDistance


def HyperbolaFromFociAndPoint(
    focus1: Matrix | Sequence[Expr],
    focus2: Matrix | Sequence[Expr],
    point: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs a hyperbola from its focus points and an incident point.

    *Formula*:
    [research/conic_from_foci_and_radius.py](../src/research/conic_from_foci_and_radius.py)
    """
    r = (PointPointDistance(focus1, point) - PointPointDistance(focus2, point)) / 2
    return ConicFromFociAndRadius(focus1, focus2, r)
