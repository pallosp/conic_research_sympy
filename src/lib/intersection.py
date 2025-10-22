from collections.abc import Sequence

from sympy import Expr, Matrix, Piecewise, nan, sqrt
from sympy.core.numbers import NaN

from lib.matrix import NonzeroCross, skew_matrix


def line_x_line(line1: Matrix, line2: Matrix) -> Matrix:
    """Computes the intersection of two lines.

    Returns an ideal point if the lines are parallel, or `[0, 0, 0]ᵀ` if they
    coincide.
    """
    return line1.cross(line2)


def conic_x_line(
    conic: Matrix,
    line: Matrix,
) -> tuple[Matrix | Sequence[Expr], Matrix | Sequence[Expr]] | NaN:
    """Intersects a conic with a line. Returns two points.

    Special cases:
     - The intersection points coincide if the line is tangent to the conic.
     - They are complex conjugates if the line doesn't intersect the conic at
       a real point.
     - Returns an unevaluated `sympy.Function` for symbolic conics.
     - Returns `None` if the conic contains the entire line.

    *Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 11.3
    """
    skew_mat = skew_matrix(line)
    m = skew_mat.T * conic * skew_mat
    a, b, c = line
    alpha = Piecewise(
        (sqrt(m[5] * m[7] - m[4] * m[8]) / a, a != 0),
        (sqrt(m[2] * m[6] - m[0] * m[8]) / b, b != 0),
        (sqrt(m[1] * m[3] - m[0] * m[4]) / c, c != 0),
    )
    if alpha == nan:
        return nan
    intersections = m + alpha * skew_mat
    points = NonzeroCross(intersections)
    if isinstance(points, (NonzeroCross, NaN)):
        return points
    return (points[0], points[1].T)
