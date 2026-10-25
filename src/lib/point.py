from sympy import I, ImmutableMatrix, Matrix, MatrixBase, S

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


def point_to_xy(point: FinitePoint) -> Matrix:
    """Computes the Euclidean coordinates of a projective point."""
    if len(point) not in (2, 3):
        raise ValueError("The point must be a 2D or 3D vector, list or tuple.")
    x, y, *z = Matrix(point).flat()
    z = z[0] if z else S.One
    return Matrix([x / z, y / z])


def point_to_vec3(point: Point) -> PointVector:
    """Computes the homogeneous coordinates of a projective point.

    Returns 3D column vectors as they are, without copying, and reshapes 3D row
    vectors to columns.
    """
    if len(point) not in (2, 3):
        raise ValueError("The point must be a 2D or 3D vector, list or tuple.")
    if len(point) == 2:
        return Matrix([point[0], point[1], 1])
    if isinstance(point, MatrixBase):
        return point if point.shape == (3, 1) else point.reshape(3, 1)
    return Matrix(point)


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
    a, b, c = line
    f = (a * x + b * y + c) / (a * a + b * b)
    return Matrix([x - a * f, y - b * f])
