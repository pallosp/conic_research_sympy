from collections.abc import Sequence

from sympy import Expr, Matrix, cos, sin

from lib.matrix import ConicMatrix
from lib.point import PointToXY


def Ellipse(
    center: Matrix | Sequence[Expr],
    r1: Expr,
    r2: Expr,
    *,
    r1_angle: Expr = None,
    r1_direction: Expr = None,
) -> Matrix:
    """Constructs an ellipse from its center, radii and the rotation angle of
    the first radius.

    *Formula*:
    [research/ellipse_from_params.py](../src/research/ellipse_from_params.py)
    """
    if (r1_angle is not None) and (r1_direction is not None):
        raise ValueError("Specify either r1_angle or r1_direction, not both.")
    center_x, center_y = PointToXY(center)
    axis_dir_x, axis_dir_y = (
        r1_direction or (1, 0) if r1_angle is None else (cos(r1_angle), sin(r1_angle))
    )
    a = -((r1 * axis_dir_y) ** 2) - (r2 * axis_dir_x) ** 2
    b = (r1**2 - r2**2) * axis_dir_x * axis_dir_y
    c = -((r1 * axis_dir_x) ** 2) - (r2 * axis_dir_y) ** 2
    d = -a * center_x - b * center_y
    e = -b * center_x - c * center_y
    f = (r1**2 * r2**2) * (axis_dir_x**2 + axis_dir_y**2) - d * center_x - e * center_y
    return ConicMatrix(a, b, c, d, e, f)


def SteinerEllipse(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
    point3: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs the Steiner circumellipse for the given points.

    The ellipse goes through the three points and is centered at the triangle's
    centroid.

    *Definition*: https://en.wikipedia.org/wiki/Steiner_ellipse<br>
    *Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)
    """
    x1, y1 = PointToXY(point1)
    x2, y2 = PointToXY(point2)
    x3, y3 = PointToXY(point3)
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
    return ConicMatrix(a, b, c, d, e, f)


def SteinerInellipse(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
    point3: Matrix | Sequence[Expr],
) -> Matrix:
    """Computes the Steiner inellipse for the given points.

    The ellipse is centered at the triangle's centroid, and is tangent to the
    triangle's sides at their midpoints.

    *Definition*: https://en.wikipedia.org/wiki/Steiner_inellipse<br>
    *Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)
    """
    x1, y1 = PointToXY(point1)
    x2, y2 = PointToXY(point2)
    x3, y3 = PointToXY(point3)
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
    return ConicMatrix(a, b, c, d, e, f)
