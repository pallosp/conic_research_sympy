from collections.abc import Sequence

from sympy import Expr, Matrix

from lib.matrix import quadratic_form, skew_matrix
from lib.point import point_to_vec3


def line_contains_point(line: Matrix, point: Matrix | Sequence[Expr]) -> bool | None:
    """Tells whether `point` is on `line`.

    Returns `None` if undecidable.
    """
    return line.dot(point_to_vec3(point)).expand().is_zero


def conic_contains_point(conic: Matrix, point: Matrix | Sequence[Expr]) -> bool | None:
    """Checks if a point lies on a conic.

    Returns `None` if undecidable.
    """
    return quadratic_form(conic, point_to_vec3(point)).expand().is_zero


def conic_contains_line(conic: Matrix, line: Matrix) -> bool | None:
    """Checks if a line lies on a conic.

    Returns `None` if undecidable.

    *Formula*:
    [research/conic_line_containment.py](../src/research/conic_line_containment.py)
    """
    skew = skew_matrix(line)
    return (skew * conic * skew).is_zero_matrix


def are_collinear(*points: Matrix) -> bool | None:
    """Tells whether n points are collinear.

    Returns `None` if undecidable.
    """
    if len(points) <= 2:
        return True
    # Treat projective points as lines.
    return are_concurrent(*(point_to_vec3(p) for p in points))


def are_concurrent(*lines: Matrix) -> bool | None:
    """Tells whether n lines are concurrent, i.e. go through the same point.

    Returns `None` if undecidable.

    Leverages the projective point-line duality, and uses the collinearity
    formula described at
    https://en.wikipedia.org/wiki/Incidence_(geometry)#Collinearity
    """
    if len(lines) <= 2:
        return True
    lines_as_matrix = Matrix.hstack(*lines)
    if len(lines) == 3:
        return lines_as_matrix.det().is_zero
    return (lines_as_matrix * lines_as_matrix.T).det().is_zero
