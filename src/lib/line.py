from typing import Sequence
from sympy import Expr, Matrix

from lib.point import PointToVec3, PointToXY


def HorizontalLine(y: Expr) -> Matrix:
    """Constructs a horizontal line with the given y-coordinate."""
    return Matrix([0, 1, -y])


def VerticalLine(x: Expr) -> Matrix:
    """Constructs a vertical line with the given x-coordinate."""
    return Matrix([-1, 0, x])


def LineBetween(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
) -> Matrix:
    """Connects two projective points with a line."""
    return PointToVec3(point1).cross(PointToVec3(point2))


def ParallelLine(
    withLine: Matrix,
    throughPoint: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs a line through a point parallel to a line."""
    x, y = PointToXY(throughPoint)
    a, b, _ = withLine
    return Matrix([a, b, -a * x - b * y])


def PerpendicularLine(
    toLine: Matrix,
    throughPoint: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs a line through a point perpendicular to a line."""
    x, y = PointToXY(throughPoint)
    a, b, _ = toLine
    return Matrix([-b, a, -b * x + a * y])


def PerpendicularBisector(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs the perpendicular bisector of two points."""
    x1, y1 = PointToXY(point1)
    x2, y2 = PointToXY(point2)
    return Matrix([x1 - x2, y1 - y2, (x2**2 - x1**2 + y2**2 - y1**2) / 2])


def AreParallel(line1: Matrix, line2: Matrix) -> bool | None:
    """Tells whether line1 and line2 are parallel.

    Returns True if they are parallel, False if not, None if undecidable.
    Considers the ideal line parallel to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * b2 - a2 * b1).expand().is_zero


def ArePerpendicular(line1: Matrix, line2: Matrix) -> bool | None:
    """Tells whether line1 and line2 are perpendicular.

    Returns True if they are perpendicular, False if not, None if undecidable.
    Considers the ideal line perpendicular to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * a2 + b1 * b2).expand().is_zero


IDEAL_LINE = Matrix([0, 0, 1])
X_AXIS = HorizontalLine(0)
Y_AXIS = VerticalLine(0)
