from collections.abc import Sequence

from sympy import Expr, Matrix

from lib.central_conic import conic_from_foci_and_radius
from lib.distance import point_point_distance


def hyperbola_from_foci_and_point(
    focus1: Matrix | Sequence[Expr],
    focus2: Matrix | Sequence[Expr],
    point: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs a hyperbola from its focus points and an incident point.

    *Formula*:
    [research/conic_from_foci_and_radius.py](../src/research/conic_from_foci_and_radius.py)
    """
    r = (point_point_distance(focus1, point) - point_point_distance(focus2, point)) / 2
    return conic_from_foci_and_radius(focus1, focus2, r)
