from collections.abc import Sequence

from sympy import Expr, Matrix, sqrt

from lib.point import PointToVec3, PointToXY


def HorizontalLine(y: Expr) -> Matrix:
    """Constructs a horizontal line with the given y-coordinate."""
    return Matrix([0, 1, -y])


def VerticalLine(x: Expr) -> Matrix:
    """Constructs a vertical line with the given x-coordinate."""
    return Matrix([-1, 0, x])


def LineBetween(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
) -> Matrix:
    """Connects two projective points with a line."""
    return PointToVec3(point1).cross(PointToVec3(point2))


def ParallelLine(
    with_line: Matrix,
    through_point: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs a line through a point parallel to a line."""
    x, y = PointToXY(through_point)
    a, b, _ = with_line
    return Matrix([a, b, -a * x - b * y])


def PerpendicularLine(
    to_line: Matrix,
    through_point: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs a line through a point perpendicular to a line."""
    x, y = PointToXY(through_point)
    a, b, _ = to_line
    return Matrix([-b, a, -b * x + a * y])


def LineThroughPoint(
    point: Matrix | Sequence[Expr],
    *,
    direction: Matrix | Sequence[Expr] = None,
    normal: Matrix | Sequence[Expr] = None,
) -> Matrix:
    """Constructs a line through a point with the given direction.

    The direction can be specified as
     - a 2D direction vector: `direction=(dx, dy)`
     - an ideal point on the line: `direction=(dx, dy, 0)`
     - a 2D normal vector: `normal=(nx, ny)`
     - an ideal point on the perpendicular line: `normal=(nx, ny, 0)`
    """
    # Exactly one of direction or normal must be specified
    if (direction is None) == (normal is None):
        raise ValueError("Specify exactly one of direction or normal.")
    x, y = PointToXY(point)
    if direction is not None:
        dx, dy, *rest = direction
        if rest not in ([], [0]):
            raise ValueError("Invalid direction vector.")
        return Matrix([-dy, dx, dy * x - dx * y])
    nx, ny, *rest = normal
    if rest not in ([], [0]):
        raise ValueError("Invalid normal vector.")
    return Matrix([nx, ny, -nx * x - ny * y])


def AngleBisector(line1: Matrix, line2: Matrix) -> Matrix:
    """Constructs the angle bisector of two lines.

    Chooses the bisector whose points substituted into the lines' equations have
    the same sign. Negate one of the lines to get the other angle bisector.

    Special cases:
     - AngleBisector(parallel finite lines, opposite direction) = center line
     - AngleBisector(coincident finite lines, same direction) = zero vector
     - AngleBisector(other parallel finite lines, same direction) = ideal line
     - AngleBisector(ideal line, ideal line) = zero vector
     - AngleBisector(ideal line, finite line) = ideal line

    *Formula*: [research/angle_bisector.py](../src/research/angle_bisector.py)
    """
    l1 = sqrt(line1[0] ** 2 + line1[1] ** 2)
    l2 = sqrt(line2[0] ** 2 + line2[1] ** 2)
    return line1 * l2 - line2 * l1


def PerpendicularBisector(
    point1: Matrix | Sequence[Expr],
    point2: Matrix | Sequence[Expr],
) -> Matrix:
    """Constructs the perpendicular bisector of two points."""
    x1, y1 = PointToXY(point1)
    x2, y2 = PointToXY(point2)
    return Matrix([x1 - x2, y1 - y2, (x2**2 - x1**2 + y2**2 - y1**2) / 2])


def AreParallel(line1: Matrix, line2: Matrix) -> bool | None:
    """Tells whether line1 and line2 are parallel.

    Returns True if they are parallel, False if not, None if undecidable.
    Considers the ideal line parallel to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * b2 - a2 * b1).expand().is_zero


def ArePerpendicular(line1: Matrix, line2: Matrix) -> bool | None:
    """Tells whether line1 and line2 are perpendicular.

    Returns True if they are perpendicular, False if not, None if undecidable.
    Considers the ideal line perpendicular to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * a2 + b1 * b2).expand().is_zero


#: The projective line at infinity.
IDEAL_LINE = Matrix([0, 0, 1])

#: The horizontal line at y=0. Substituting a point at y>0 to the line's
#: equation will evaluate to a positive value.
X_AXIS = HorizontalLine(0)

#: The vertical line at x=0. Substituting a point at x<0 to the line's equation
#: will evaluate to a positive value.
Y_AXIS = VerticalLine(0)
