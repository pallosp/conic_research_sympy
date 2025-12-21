from collections.abc import Sequence

from sympy import Expr, Matrix, sqrt

from lib.matrix import conic_matrix
from lib.point import ORIGIN, point_to_xy


def circle(center: Matrix | Sequence[Expr], radius: Expr) -> Matrix:
    """Creates a circle from its center and radius."""
    x, y = point_to_xy(center)
    return conic_matrix(-1, 0, -1, x, y, radius * radius - x * x - y * y)


def circle_radius(circle: Matrix) -> Expr:
    """Computes the radius of a circle conic.

    The result is not specified if the conic matrix is not a circle.
    The computation is based on
    [research/director_circle.py](../src/research/director_circle.py).
    """
    a, b, c = circle[0], circle[1], circle[4]
    return sqrt(-circle.det() * (a + c) / 2) / (a * c - b * b)


def director_circle(conic: Matrix) -> Matrix:
    """Computes the director circle of a conic. It's also called orthoptic
    circle or Fermat–Apollonius circle.

    *Definition*: https://en.wikipedia.org/wiki/Director_circle<br>
    *Formula*: [research/director_circle.py](../src/research/director_circle.py)
    """
    a, _, _, _, c, _, d, e, f = conic.adjugate()
    return Matrix(
        [
            [-1, 0, d / f],
            [0, -1, e / f],
            [d / f, e / f, -(a + c) / f],
        ],
    )


#: The circle at the origin with radius 1.
UNIT_CIRCLE: Matrix = circle(ORIGIN, 1)

#: The circle at the origin with radius 𝑖.
IMAGINARY_UNIT_CIRCLE: Matrix = Matrix.eye(3)
