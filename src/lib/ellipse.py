from collections.abc import Sequence

from sympy import Expr, Matrix, cos, sin

from lib.central_conic import ConicFromFociAndRadius
from lib.distance import PointPointDistance
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
    """Constructs an ellipse from its center, radii, and the either the
    direction vector of the first radius or its angle to horizontal.

    *Formula*:
    [research/ellipse_from_params.py](../src/research/ellipse_from_params.py)
    """
    if (r1_angle is not None) and (r1_direction is not None):
        raise ValueError("Specify either r1_angle or r1_direction, not both.")

    # center
    cx, cy = PointToXY(center)

    # major axis direction
    dx, dy = 1, 0
    if r1_angle is not None:
        dx, dy = cos(r1_angle), sin(r1_angle)
    elif r1_direction is not None:
        dx, dy, *rest = r1_direction
        if rest not in ([], [0]):
            raise ValueError("r1_direction must be a 2d vector or an ideal point")

    a = -((r1 * dy) ** 2) - (r2 * dx) ** 2
    b = (r1**2 - r2**2) * dx * dy
    c = -((r1 * dx) ** 2) - (r2 * dy) ** 2
    d = -a * cx - b * cy
    e = -b * cx - c * cy
    f = (r1**2 * r2**2) * (dx**2 + dy**2) - d * cx - e * cy

    return ConicMatrix(a, b, c, d, e, f)


def EllipseFromFociAndPoint(
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
    [research/ellipse_from_foci_and_point.py](../src/research/ellipse_from_foci_and_point.py)
    """
    r = (PointPointDistance(focus1, point) + PointPointDistance(focus2, point)) / 2
    return ConicFromFociAndRadius(focus1, focus2, r)


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
