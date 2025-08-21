from typing import Sequence
from sympy import Expr, Matrix, sqrt

from lib.matrix import ConicMatrix
from lib.point import ORIGIN, PointToXY


def Circle(center: Matrix | Sequence[Expr], radius: Expr) -> Matrix:
    """Creates a circle from its center and radius."""
    x, y = PointToXY(center)
    return ConicMatrix(-1, 0, -1, x, y, radius * radius - x * x - y * y)


def CircleRadius(circle: Matrix) -> Expr:
    """Computes the radius of a circle conic.

    The result is not specified if the conic matrix is not a circle.
    The computation is based on `research/director_circle.py`.
    """
    a, b, c = circle[0], circle[1], circle[4]
    return sqrt(-circle.det() * (a + c) / 2) / (a * c - b * b)


def DirectorCircle(conic: Matrix) -> Matrix:
    """Computes the director circle of a conic. It's also called orthoptic
    circle or Fermat–Apollonius circle.

    *Definition*: https://en.wikipedia.org/wiki/Director_circle<br>
    *Formula*: `research/director_circle.py`
    """
    a, _, _, _, c, _, d, e, f = conic.adjugate()
    return Matrix(
        [
            [-1, 0, d / f],
            [0, -1, e / f],
            [d / f, e / f, -(a + c) / f],
        ]
    )


#: The circle at the origin with radius 1.
UNIT_CIRCLE = Circle(ORIGIN, 1)

#: The circle at the origin with radius 𝑖.
COMPLEX_UNIT_CIRCLE = Matrix.eye(3)
