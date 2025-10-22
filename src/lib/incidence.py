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

    Algorithm:
     - n=3: https://en.wikipedia.org/wiki/Incidence_(geometry)#Collinearity
     - n>3: https://en.wikipedia.org/wiki/Gram_matrix
    """
    if len(points) <= 2:
        return True
    points_as_matrix = Matrix.hstack(*(point_to_vec3(p) for p in points))
    if len(points) == 3:
        return points_as_matrix.det().is_zero
    return (points_as_matrix * points_as_matrix.T).det().is_zero
