from typing import Sequence
from sympy import Expr, Matrix, sqrt

from lib.point import PointToXY, PointToVec3


def PointPointDistance(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
) -> Expr:
    """Computes the signed distance between two points.

    Special cases:
     - the distance between Euclidean and ideal points is infinity
     - the distance between two ideal points is `nan`
    """
    x1, y1, z1 = PointToVec3(point1)
    x2, y2, z2 = PointToVec3(point2)
    return sqrt((x2 * z1 - x1 * z2) ** 2 + (y2 * z1 - y1 * z2) ** 2) / (z1 * z2)


def PointLineDistance(point: Matrix | Sequence[Expr], line: Matrix) -> Expr:
    """Computes the signed distance between a point and a line.

    Special cases:
     - the distance between Euclidean points and an ideal lines is infinity
     - the distance between ideal points and an any projective lines is `nan`
    """
    x, y = PointToXY(point)
    a, b, c = line
    return (a * x + b * y + c) / sqrt(a * a + b * b)
