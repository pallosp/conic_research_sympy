from sympy import Matrix

from lib.matrix import conic_matrix
from lib.point import point_to_vec3, point_to_xy
from lib.type_aliases import FinitePoint, Point, PointVector


def circumcircle(
    point1: Point,
    point2: Point,
    point3: Point,
) -> Matrix:
    """Computes the circle through three points.

    If the points are in counterclockwise order, the sign of the matrix matches
    `circle()`: the conic's value at the center is non-negative. In clockwise
    order, the matrix is negated.

    Special cases:
     - If the points are collinear, the result is the union of their common
       line and the ideal line.
     - If exactly one point is ideal, the result is the union of the line
       through the two finite points and the ideal line.
     - If two points coincide, or at least two points are ideal, the result is
       a zero matrix.

    *Definition*: <https://en.wikipedia.org/wiki/Circumcircle><br>
    *Formula*:
    [research/construction/circumcircle.py](../src/research/construction/circumcircle.py)
    """
    x1, y1, z1 = point_to_vec3(point1)
    x2, y2, z2 = point_to_vec3(point2)
    x3, y3, z3 = point_to_vec3(point3)
    s_col = [x1**2 + y1**2, x2**2 + y2**2, x3**2 + y3**2]
    x_col = [x1, x2, x3]
    y_col = [y1, y2, y3]
    z_col = [z1, z2, z3]
    xz_col = [x1 * z1, x2 * z2, x3 * z3]
    yz_col = [y1 * z1, y2 * z2, y3 * z3]
    zz_col = [z1**2, z2**2, z3**2]
    a = -2 * z1 * z2 * z3 * Matrix([x_col, y_col, z_col]).det()
    d = Matrix([s_col, yz_col, zz_col]).det()
    e = -Matrix([s_col, xz_col, zz_col]).det()
    f = 2 * Matrix([s_col, xz_col, yz_col]).det()
    return conic_matrix(a, 0, a, d, e, f)


def circumcenter(
    point1: Point,
    point2: Point,
    point3: Point,
) -> PointVector:
    """Computes the center of the circle through three points.

    Returns the center as a homogeneous 3D column vector.

    Special cases:
     - If the points are collinear, the result is the ideal point
       perpendicular to their common line.
     - If exactly one point is ideal, the result is the ideal point
       perpendicular to the line through the two finite points.
     - If two points coincide, or at least two points are ideal, the result is
       a zero vector.

    *Definition*: <https://en.wikipedia.org/wiki/Circumcircle><br>
    *Formula*:
    [research/construction/circumcircle.py](../src/research/construction/circumcircle.py)
    """
    x1, y1, z1 = point_to_vec3(point1)
    x2, y2, z2 = point_to_vec3(point2)
    x3, y3, z3 = point_to_vec3(point3)
    s_col = [x1**2 + y1**2, x2**2 + y2**2, x3**2 + y3**2]
    x_col = [x1, x2, x3]
    y_col = [y1, y2, y3]
    z_col = [z1, z2, z3]
    xz_col = [x1 * z1, x2 * z2, x3 * z3]
    yz_col = [y1 * z1, y2 * z2, y3 * z3]
    zz_col = [z1**2, z2**2, z3**2]
    x = Matrix([s_col, yz_col, zz_col]).det()
    y = -Matrix([s_col, xz_col, zz_col]).det()
    z = 2 * z1 * z2 * z3 * Matrix([x_col, y_col, z_col]).det()
    return Matrix([x, y, z])


def steiner_ellipse(
    point1: FinitePoint,
    point2: FinitePoint,
    point3: FinitePoint,
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
    point1: FinitePoint,
    point2: FinitePoint,
    point3: FinitePoint,
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
