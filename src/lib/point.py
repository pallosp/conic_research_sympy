from typing import Sequence, Tuple
from sympy import Expr, Matrix, sympify

ORIGIN = Matrix([0, 0, 1])


def IdealPoint(x: Expr, y: Expr) -> Matrix:
    """Creates an ideal point at the given direction."""
    return Matrix([x, y, 0])


def IdealPointOnLine(line: Matrix) -> Matrix:
    """Returns the coordinates of the ideal point on the line.

    The first two coordinates specify the line's direction. The third one is
    always zero.

    If the line is the ideal line, returns a zero vector.
    """
    return Matrix([line[1], -line[0], 0])


def PointToXY(point: Matrix | Sequence[Expr]) -> Tuple[Expr, Expr]:
    """Computes the Euclidean coordinates of a projective point."""
    assert len(point) in (2, 3)
    if len(point) == 2:
        return sympify(point)
    x, y, z = sympify(point)
    return (x / z, y / z)


def PointToVec3(point: Matrix | Sequence[Expr]) -> Matrix:
    """Computes the homogeneous coordinates of a projective point."""
    assert len(point) in (2, 3)
    if len(point) == 2:
        return Matrix([point[0], point[1], 1])
    if point is Matrix:
        return point
    return Matrix(point)


def Centroid(*points: Sequence[Expr]) -> Tuple[Expr, Expr]:
    """Computes the centroid of a set of points."""
    n = len(points)
    assert n > 0
    cx, cy = sympify((0, 0))
    for p in points:
        x, y = PointToXY(p)
        cx += x
        cy += y
    return (cx / n, cy / n)
