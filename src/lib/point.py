from collections.abc import Sequence

from sympy import Expr, Matrix, S, sympify

#: The point at (0, 0)
ORIGIN: Matrix = Matrix([0, 0, 1])


def ideal_point(x: Expr, y: Expr) -> Matrix:
    """Creates an ideal point at the given direction."""
    return Matrix([x, y, 0])


def ideal_point_on_line(line: Matrix) -> Matrix:
    """Returns the coordinates of the ideal point on the line.

    The first two coordinates specify the line's direction. The third one is
    always zero.

    If the line is the ideal line, returns `[0, 0, 0]ᵀ`.
    """
    return Matrix([line[1], -line[0], 0])


def point_to_xy(point: Matrix | Sequence[Expr]) -> tuple[Expr, Expr]:
    """Computes the Euclidean coordinates of a projective point."""
    if len(point) not in (2, 3):
        raise ValueError("The point must be a 2D or 3D vector, list or tuple.")
    if len(point) == 2:
        return sympify(point)
    x, y, z = sympify(point)
    return (x / z, y / z)


def point_to_vec3(point: Matrix | Sequence[Expr]) -> Matrix:
    """Computes the homogeneous coordinates of a projective point."""
    if len(point) not in (2, 3):
        raise ValueError("The point must be a 2D or 3D vector, list or tuple.")
    if len(point) == 2:
        return Matrix([point[0], point[1], 1])
    if point is Matrix:
        return point
    return Matrix(point)


def centroid(*points: Sequence[Expr]) -> Matrix:
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
    point: Matrix | Sequence[Expr],
    line: Matrix,
) -> Matrix:
    """Computes the foot of the perpendicular through a point to a line.

    Returns the foot point's coordinates as a 2D column vector.
    Degenerates to `[nan, nan]ᵀ` when `point`, `line` or both are infinite.

    *Formula*: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    """
    x, y = point_to_xy(point)
    a, b, c = line
    f = (a * x + b * y + c) / (a * a + b * b)
    return Matrix([x - a * f, y - b * f])
