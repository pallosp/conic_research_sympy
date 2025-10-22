from collections.abc import Sequence

from sympy import Expr, Matrix, sqrt

from lib.line import are_parallel
from lib.point import point_to_vec3, point_to_xy


def point_point_distance(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
) -> Expr:
    """Computes the signed distance between two points.

    Special cases:
     - the distance between finite and ideal points is infinity
     - the distance between two ideal points is `nan`
    """
    x1, y1, z1 = point_to_vec3(point1)
    x2, y2, z2 = point_to_vec3(point2)
    return sqrt((x2 * z1 - x1 * z2) ** 2 + (y2 * z1 - y1 * z2) ** 2) / (z1 * z2)


def point_line_distance(point: Matrix | Sequence[Expr], line: Matrix) -> Expr:
    """Computes the signed distance between a point and a line.

    Special cases:
     - the distance between finite points and the ideal line is infinity
     - the distance between ideal points and any projective lines is `nan`
    """
    x, y = point_to_xy(point)
    a, b, c = line
    return (a * x + b * y + c) / sqrt(a * a + b * b)


def parallel_line_distance(line1: Matrix, line2: Matrix) -> Expr:
    """Computes the signed distance between two parallel lines.

    The return value is positive if the two lines have the same direction,
    negative otherwise.

    Special cases:
    - Returns infinity if one of the lines is the ideal line.
    - Returns `nan` if both lines are ideal lines.
    - Raises a `ValueError` if the lines are provably not parallel.
    - Returns an unspecified value if the lines cross at a finite point, but
      Sympy cannot prove this fact.
    """
    if are_parallel(line1, line2) is False:
        raise ValueError("The lines must be parallel")
    a1, b1, c1 = line1
    a2, b2, c2 = line2
    return sqrt((a1 * c2 - a2 * c1) ** 2 + (b1 * c2 - b2 * c1) ** 2) / (
        a1 * a2 + b1 * b2
    )
