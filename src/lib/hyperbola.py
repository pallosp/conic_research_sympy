from collections.abc import Sequence

from sympy import Expr, Matrix, acos, sqrt

from lib.central_conic import (
    conic_from_foci_and_radius,
)
from lib.conic_classes import ConicNormFactor
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


def asymptote_focal_axis_angle(hyperbola: Matrix) -> Expr:
    """Computes the angle between the axis and an asymptote of a hyperbola.

    Return values for other conics:
    - *non-degenarate conics*
      - parabolas: 0
      - real and imaginary circles: infinity
      - other real ellipses: imaginary angle (`acos(x)`; `x > 1`)
      - non-circular complex ellipses: complex angle
    - *point conics*
      - zero-radius circles: infinity
      - non-circular finite point conics: imaginary angle
      - ideal point conics: 0
    - *line pairs*
       - crossing lines: positive angle; α or π/2-α depending on the lines'
         directions
       - parallel lines in the same direction: π/2
       - parallel lines in the opposite direction: 0
       - conics containis the ideal line: `nan`

    *Research*:
    [research/asymptote_angle.py](../src/research/asymptote_angle.py)<br>
    *Formula*: `acos(1 / eccentricity)`
    """
    a, _, _, b, c, _, _, _, _ = hyperbola
    s = sqrt(((a - c) ** 2 + 4 * b**2).factor())
    norm_sign = ConicNormFactor(hyperbola)
    return acos(sqrt((s - norm_sign * (a + c)) / (2 * s)))
