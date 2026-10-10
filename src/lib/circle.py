from collections.abc import Sequence

from sympy import Expr, Matrix, sqrt

from lib.matrix import conic_matrix
from lib.point import ORIGIN, point_to_vec3, point_to_xy


def circle(
    center: Matrix | Sequence[Expr],
    *,
    r: Expr | None = None,
    point: Matrix | Sequence[Expr] | None = None,
) -> Matrix:
    """Creates a circle from its center and either its radius or a point on it.

    Exactly one of `r` and `point` must be specified.

    The conic's value at the center is non-negative (`r²` for a finite
    radius `r`), so the sign of the matrix is the same for both constructions.

    If `center` is an ideal point, the result is a matrix with `nan` and/or
    `zoo` entries. If `point` is an ideal point and the center is finite, the
    result is the double ideal line.

    *Formula*:
    [research/construction/circle.py](../src/research/construction/circle.py)
    """
    if (r is None) == (point is None):
        raise ValueError("Exactly one of r and point must be specified.")
    cx, cy = point_to_xy(center)
    if r is None:
        x, y, z = point_to_vec3(point)
        return conic_matrix(
            -(z**2),
            0,
            -(z**2),
            cx * z**2,
            cy * z**2,
            x**2 + y**2 - 2 * z * (cx * x + cy * y),
        )
    return conic_matrix(-1, 0, -1, cx, cy, r * r - cx * cx - cy * cy)


def circle_radius(conic: Matrix) -> Expr:
    """Computes the radius of a circular or a finite point conic.

    Return value by conic type:

    - *Real circles*: the positive radius.
    - *Imaginary circles*: an imaginary number.
    - *Finite point conics*: 0, including non-circular ones.
    - *Other conics*: not meaningful.

    *Formula*:
    [research/conic_properties/circle_radius.py](../src/research/conic_properties/circle_radius.py)
    """
    a = conic[0]
    return sqrt(-a * conic.det()) / a**2


def director_circle(conic: Matrix) -> Matrix:
    """Computes the director circle of a conic.

    It's also called orthoptic circle or Fermat–Apollonius circle.

    *Definition*: <https://en.wikipedia.org/wiki/Director_circle><br>
    *Formula*:
    [research/construction/director_circle.py](../src/research/construction/director_circle.py)
    """
    a, _, _, _, c, _, d, e, f = conic.adjugate()
    return Matrix(
        [
            [-1, 0, d / f],
            [0, -1, e / f],
            [d / f, e / f, -(a + c) / f],
        ],
    )


#: The circle at the origin with radius 1.
UNIT_CIRCLE: Matrix = circle(ORIGIN, r=1)

#: The circle at the origin with radius 𝑖.
IMAGINARY_UNIT_CIRCLE: Matrix = Matrix.eye(3)
