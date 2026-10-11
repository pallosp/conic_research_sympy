from collections.abc import Sequence

from sympy import Expr, Matrix, cos, sin

from lib.central_conic import conic_from_foci_and_radius
from lib.distance import point_point_distance
from lib.matrix import conic_matrix
from lib.point import point_to_xy


def ellipse(
    center: Matrix | Sequence[Expr],
    r1: Expr,
    r2: Expr,
    *,
    r1_angle: Expr = None,
    r1_direction: Expr = None,
) -> Matrix:
    """Constructs an ellipse from its center, radii, and the either the
    direction vector of the first radius or its angle to horizontal.

    *Formula*:
    [research/construction/ellipse_from_params.py](../src/research/construction/ellipse_from_params.py)
    """
    if (r1_angle is not None) and (r1_direction is not None):
        raise ValueError("Specify either r1_angle or r1_direction, not both.")

    # center
    cx, cy = point_to_xy(center)

    # major axis direction
    dx, dy = 1, 0
    if r1_angle is not None:
        dx, dy = cos(r1_angle), sin(r1_angle)
    elif r1_direction is not None:
        dx, dy, *rest = r1_direction
        if rest not in ([], [0]):
            raise ValueError("r1_direction must be a 2D vector or an ideal point")

    a = -((r1 * dy) ** 2 + (r2 * dx) ** 2)
    b = (r1**2 - r2**2) * dx * dy
    c = -((r1 * dx) ** 2 + (r2 * dy) ** 2)
    d = -a * cx - b * cy
    e = -b * cx - c * cy
    f = r1**2 * r2**2 * (dx**2 + dy**2) - d * cx - e * cy

    return conic_matrix(a, b, c, d, e, f)


def ellipse_from_foci_and_point(
    focus1: Matrix | Sequence[Expr],
    focus2: Matrix | Sequence[Expr],
    point: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs an ellipse from its focus points and an incident point.

    Special cases:
     - If all three points conicide, returns a zero matrix.
     - If the foci coincide but the point doesn't the result is a circle.
     - If the points are collinear but the foci don't coincide, returns a
       coincident line conic.
     - If any of the points are ideal points, returns a matrix that contains
       `nan` elements.

    *Formula*:
    [research/construction/ellipse_from_foci_and_point.py](../src/research/construction/ellipse_from_foci_and_point.py)
    """
    r = (point_point_distance(focus1, point) + point_point_distance(focus2, point)) / 2
    return conic_from_foci_and_radius(focus1, focus2, r)
