from collections.abc import Sequence

from sympy import ImmutableMatrix, Matrix, MatrixBase, sign, sqrt

from lib.point import direction_to_xy, ideal_point, point_to_vec3, point_to_xy
from lib.type_aliases import FinitePoint, Line, Point, PointVector, Scalar


def horizontal_line(y: Scalar) -> Line:
    """Constructs a horizontal line with the given y-coordinate."""
    return Matrix([0, 1, -y])


def vertical_line(x: Scalar) -> Line:
    """Constructs a vertical line with the given x-coordinate."""
    return Matrix([-1, 0, x])


def line_between(
    point1: Point,
    point2: Point,
) -> Line:
    """Connects two projective points with a line."""
    return point_to_vec3(point1).cross(point_to_vec3(point2))


def parallel_line(
    with_line: Line,
    through_point: FinitePoint,
) -> Line:
    """Constructs a line through a point parallel to a line."""
    x, y = point_to_xy(through_point)
    a, b, _ = with_line
    return Matrix([a, b, -a * x - b * y])


def perpendicular_line(
    to_line: Line,
    through_point: FinitePoint,
) -> Line:
    """Constructs a line through a point perpendicular to a line."""
    x, y = point_to_xy(through_point)
    a, b, _ = to_line
    return Matrix([-b, a, -b * x + a * y])


def line_through_point(
    point: FinitePoint,
    *,
    direction: MatrixBase | Sequence[Scalar] | None = None,
    normal: MatrixBase | Sequence[Scalar] | None = None,
) -> Line:
    """Constructs a line through a point with the given direction.

    The direction can be specified as
     - a 2D direction vector: `direction=(dx, dy)`
     - an ideal point on the line: `direction=(dx, dy, 0)`
     - a 2D normal vector: `normal=(nx, ny)`
     - an ideal point on the perpendicular line: `normal=(nx, ny, 0)`
    """
    x, y = point_to_xy(point)
    if direction is not None and normal is None:
        dx, dy = direction_to_xy(direction)
        return Matrix([-dy, dx, dy * x - dx * y])
    if normal is not None and direction is None:
        nx, ny = direction_to_xy(normal)
        return Matrix([nx, ny, -nx * x - ny * y])
    raise ValueError("Specify exactly one of direction or normal.")


def angle_bisector(line1: Line, line2: Line) -> Line:
    """Constructs the angle bisector of two lines.

    Chooses the bisector whose points substituted into the lines' equations have
    the same sign. Negate one of the lines to get the other angle bisector.

    Special cases:
     - AngleBisector(parallel finite lines, opposite direction) = center line
     - AngleBisector(coincident finite lines, same direction) = `[0, 0, 0]ᵀ`
     - AngleBisector(other parallel finite lines, same direction) = ideal line
     - AngleBisector(ideal line, ideal line) = `[0, 0, 0]ᵀ`
     - AngleBisector(ideal line, finite line) = ideal line

    *Formula*:
    [research/construction/angle_bisector.py](../src/research/construction/angle_bisector.py)
    """
    l1 = sqrt(line1[0] ** 2 + line1[1] ** 2)
    l2 = sqrt(line2[0] ** 2 + line2[1] ** 2)
    return line1 * l2 - line2 * l1


def perpendicular_bisector(
    point1: FinitePoint,
    point2: FinitePoint,
) -> Line:
    """Constructs the perpendicular bisector of two points."""
    x1, y1 = point_to_xy(point1)
    x2, y2 = point_to_xy(point2)
    return Matrix([x1 - x2, y1 - y2, (x2**2 - x1**2 + y2**2 - y1**2) / 2])


def line_normal(
    line: Line,
    *,
    toward: Point | None = None,
) -> PointVector:
    """Returns the normal vector of a line, represented as an ideal point.

    By default, the normal is chosen such that `line·normal > 0`.
    If a `toward` point is provided, the returned vector points toward it in the
    projective sense: `line·normal` and `line·toward` have the same sign.
    For an ordinary Euclidean point (given by its `(x, y)` coordinates), this
    corresponds to the normal vector pointing toward the semiplane containing that
    point.

    Returns `[0, 0, 0]ᵀ` if the `toward` point lies on the line.
    """
    normal = ideal_point(line[0], line[1])
    if toward is None:
        return normal
    return normal * sign(line.dot(point_to_vec3(toward)).factor())


def are_parallel(line1: Line, line2: Line) -> bool | None:
    """Tells whether line1 and line2 are parallel.

    Returns `True` if they are parallel, `False` if not, `None` if undecidable.
    Considers the ideal line parallel to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * b2 - a2 * b1).expand().is_zero


def are_perpendicular(line1: Line, line2: Line) -> bool | None:
    """Tells whether line1 and line2 are perpendicular.

    Returns `True` if they are perpendicular, `False` if not, `None` if
    undecidable. Considers the ideal line perpendicular to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * a2 + b1 * b2).expand().is_zero


#: The zero vector `(0 : 0 : 0)`. It doesn't represent any line. Functions
#: return it when the resulting line is undefined.
ZERO_LINE: ImmutableMatrix = ImmutableMatrix([0, 0, 0])

#: The projective line at infinity.
IDEAL_LINE: ImmutableMatrix = ImmutableMatrix([0, 0, 1])

#: The horizontal line at y=0. Substituting a point at y>0 to the line's
#: equation will evaluate to a positive value.
X_AXIS: ImmutableMatrix = ImmutableMatrix(horizontal_line(0))

#: The vertical line at x=0. Substituting a point at x<0 to the line's equation
#: will evaluate to a positive value.
Y_AXIS: ImmutableMatrix = ImmutableMatrix(vertical_line(0))
