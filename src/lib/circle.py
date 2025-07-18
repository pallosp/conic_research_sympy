from sympy import Matrix, sqrt
from lib.matrix import ConicMatrix
from lib.point import ORIGIN, PointToXY


def Circle(center, r):
    x, y = PointToXY(center)
    return ConicMatrix(-1, 0, -1, x, y, r * r - x * x - y * y)


def CircleRadius(circle):
    """Computes the radius of a circle conic.

    The result is not specified if the conic matrix is not a circle.

    The computation is based on director_circle.py.
    """
    a, b, c = circle[0], circle[1], circle[4]
    return sqrt(-circle.det() * (a + c) / 2) / (a * c - b * b)


def DirectorCircle(conic):
    a, _, _, _, c, _, d, e, f = conic.adjugate()
    return Matrix(
        [
            [-1, 0, d / f],
            [0, -1, e / f],
            [d / f, e / f, -(a + c) / f],
        ]
    )


UNIT_CIRCLE = Circle(ORIGIN, 1)
