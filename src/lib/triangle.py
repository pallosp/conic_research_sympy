from collections.abc import Sequence

from sympy import Expr, Matrix

from lib.matrix import conic_matrix
from lib.point import point_to_xy


def steiner_ellipse(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
    point3: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs the Steiner circumellipse for the given points.

    The ellipse goes through the three points and is centered at the triangle's
    centroid.

    *Definition*: <https://en.wikipedia.org/wiki/Steiner_ellipse><br>
    *Formula*:
    [research/construction/steiner_ellipse.py](../src/research/construction/steiner_ellipse.py)
    """
    x1, y1 = point_to_xy(point1)
    x2, y2 = point_to_xy(point2)
    x3, y3 = point_to_xy(point3)
    x_row = [x1, x2, x3]
    y_row = [y1, y2, y3]
    dx_row = [x2 - x3, x3 - x1, x1 - x2]
    dy_row = [y2 - y3, y3 - y1, y1 - y2]
    a = Matrix([dy_row, y_row, [1, 1, 1]]).det()
    b = Matrix([x_row, dy_row, [1, 1, 1]]).det()
    c = Matrix([dx_row, x_row, [1, 1, 1]]).det()
    d = Matrix([y_row, dy_row, x_row]).det()
    e = Matrix([x_row, dx_row, y_row]).det()
    f = (
        Matrix(
            [
                x_row,
                [x1 * y2 - x2 * y1, x2 * y3 - x3 * y2, x3 * y1 - x1 * y3],
                y_row,
            ],
        ).det()
        * 2
    )
    return conic_matrix(a, b, c, d, e, f)


def steiner_inellipse(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
    point3: Matrix | Sequence[Expr],
) -> Matrix:
    """Computes the Steiner inellipse for the given points.

    The ellipse is centered at the triangle's centroid, and is tangent to the
    triangle's sides at their midpoints.

    *Definition*: <https://en.wikipedia.org/wiki/Steiner_inellipse><br>
    *Formula*:
    [research/construction/steiner_ellipse.py](../src/research/construction/steiner_ellipse.py)
    """
    x1, y1 = point_to_xy(point1)
    x2, y2 = point_to_xy(point2)
    x3, y3 = point_to_xy(point3)
    x_row = [x1, x2, x3]
    y_row = [y1, y2, y3]
    dx_row = [x2 - x3, x3 - x1, x1 - x2]
    dy_row = [y2 - y3, y3 - y1, y1 - y2]
    a = Matrix([dy_row, y_row, [1, 1, 1]]).det()
    b = Matrix([x_row, dy_row, [1, 1, 1]]).det()
    c = Matrix([dx_row, x_row, [1, 1, 1]]).det()
    d = Matrix([y_row, dy_row, x_row]).det()
    e = Matrix([x_row, dx_row, y_row]).det()
    f = (
        Matrix(
            [
                [x2 + x3, x3 + x1, x1 + x2],
                [x2 * y3 - x3 * y2, x3 * y1 - x1 * y3, x1 * y2 - x2 * y1],
                [y2 + y3, y3 + y1, y1 + y2],
            ],
        ).det()
        / 2
    )
    return conic_matrix(a, b, c, d, e, f)
