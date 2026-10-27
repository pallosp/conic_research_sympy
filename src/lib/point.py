from collections.abc import Sequence
from typing import overload

from sympy import Expr, I, ImmutableMatrix, Matrix, MatrixBase, MatrixExpr, S, sympify

from lib.type_aliases import FinitePoint, Line, Point, PointVector, Scalar

#: The zero vector `(0 : 0 : 0)`. It doesn't represent any point. Functions
#: return it when the resulting point is undefined.
ZERO_POINT: ImmutableMatrix = ImmutableMatrix([0, 0, 0])

#: The point at (0, 0)
ORIGIN: ImmutableMatrix = ImmutableMatrix([0, 0, 1])

#: The circular point at infinity `(1 : i : 0)`, also known as a cyclic or
#: isotropic point. Every circle passes through it.
#:
#: https://en.wikipedia.org/wiki/Circular_points_at_infinity
CIRCULAR_POINT_I: ImmutableMatrix = ImmutableMatrix([1, I, 0])

#: The circular point at infinity `(1 : -i : 0)`, also known as a cyclic or
#: isotropic point. It is the complex conjugate of `CIRCULAR_POINT_I`.
#: Every circle passes through it.
CIRCULAR_POINT_J: ImmutableMatrix = ImmutableMatrix([1, -I, 0])


def ideal_point(x: Scalar, y: Scalar) -> PointVector:
    """Creates an ideal point at the given direction."""
    return Matrix([x, y, 0])


def ideal_point_on_line(line: Line) -> PointVector:
    """Returns the coordinates of the ideal point on the line.

    The first two coordinates specify the line's direction. The third one is
    always zero.

    If the line is the ideal line, returns `[0, 0, 0]ᵀ`.
    """
    return ideal_point(line[1], -line[0])


def point_to_xy(point: FinitePoint) -> tuple[Expr, Expr]:
    """Computes the Euclidean coordinates of a projective point as an `(x, y)` tuple."""
    if len(point) not in (2, 3):
        raise ValueError("The point must be a 2D or 3D vector, list or tuple.")
    x, y, *z = Matrix(point).flat()
    z = z[0] if z else S.One
    return x / z, y / z


def direction_to_xy(direction: MatrixBase | Sequence[Scalar]) -> tuple[Expr, Expr]:
    """Returns the components of a 2D direction vector or an ideal point."""
    if isinstance(direction, MatrixBase):
        dx, dy, *rest = direction.flat()
    else:
        dx, dy, *rest = map(sympify, direction)
    if rest not in ([], [0]):
        raise ValueError("The direction must be a 2D vector or an ideal point.")
    return dx, dy


@overload
def point_to_vec3(point: Point) -> PointVector: ...
@overload
def point_to_vec3(point: MatrixExpr) -> PointVector | MatrixExpr: ...
def point_to_vec3(point: Point | MatrixExpr) -> PointVector | MatrixExpr:
    """Computes the homogeneous coordinates of a projective point.

    Returns 3D column vectors as they are, without copying, and transposes 3D
    row vectors to columns. Accepts unevaluated `MatrixExpr` vectors too.
    """
    if isinstance(point, (tuple, list)):
        if len(point) == 2:
            return Matrix([point[0], point[1], 1])
        if len(point) == 3:
            return Matrix(point)
    elif point.shape in ((2, 1), (1, 2)):
        return Matrix([point[0], point[1], 1])
    elif point.shape == (3, 1):
        return point
    elif point.shape == (1, 3):
        return point.T
    raise ValueError("The point must be a 2D or 3D vector, list or tuple.")


def centroid(*points: FinitePoint) -> Matrix:
    """Computes the centroid of a set of points.

    Returns the point's coordinates as a 2D column vector.
    """
    n = len(points)
    if n == 0:
        raise ValueError("At least one point must be provided.")
    cx, cy = S.Zero, S.Zero
    for p in points:
        x, y = point_to_xy(p)
        cx += x
        cy += y
    return Matrix([cx / n, cy / n])


def perpendicular_foot(
    point: FinitePoint,
    line: Line,
) -> Matrix:
    """Computes the foot of the perpendicular through a point to a line.

    Returns the foot point's coordinates as a 2D column vector.
    Degenerates to `[nan, nan]ᵀ` when `point`, `line` or both are infinite.

    *Formula*: <https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line>
    """
    x, y = point_to_xy(point)
    a, b, c = line.flat()
    f = (a * x + b * y + c) / (a * a + b * b)
    return Matrix([x - a * f, y - b * f])
