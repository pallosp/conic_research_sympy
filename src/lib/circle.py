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
    x, y = point_to_xy(center)
    if r is None:
        px, py, pz = point_to_vec3(point)
        return conic_matrix(
            -(pz**2),
            0,
            -(pz**2),
            x * pz**2,
            y * pz**2,
            px**2 + py**2 - 2 * pz * (x * px + y * py),
        )
    return conic_matrix(-1, 0, -1, x, y, r * r - x * x - y * y)


def circle_radius(circle: Matrix) -> Expr:
    """Computes the radius of a circle conic.

    The result is not specified if the conic matrix is not a circle.
    The computation is based on
    [research/construction/director_circle.py](../src/research/construction/director_circle.py).
    """
    a, b, c = circle[0], circle[1], circle[4]
    return sqrt(-circle.det() * (a + c) / 2) / (a * c - b * b)


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
