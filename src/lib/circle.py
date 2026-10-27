from sympy import Expr, ImmutableMatrix, Matrix, sqrt

from lib.matrix import conic_coefficients, conic_matrix
from lib.point import ORIGIN, point_to_vec3, point_to_xy
from lib.type_aliases import Conic, FinitePoint, Point, Scalar


def circle(
    center: FinitePoint,
    *,
    r: Scalar | None = None,
    point: Point | None = None,
) -> Conic:
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
    cx, cy = point_to_xy(center)
    if r is not None and point is None:
        return conic_matrix(-1, 0, -1, cx, cy, r * r - cx * cx - cy * cy)
    if point is not None and r is None:
        x, y, z = point_to_vec3(point).flat()
        return conic_matrix(
            -(z**2),
            0,
            -(z**2),
            cx * z**2,
            cy * z**2,
            x**2 + y**2 - 2 * z * (cx * x + cy * y),
        )
    raise ValueError("Exactly one of r and point must be specified.")


def circle_radius(conic: Conic) -> Expr:
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


def director_circle(conic: Conic) -> Conic:
    """Computes the director circle of a conic.

    It's also called orthoptic circle or Fermat–Apollonius circle. It's the
    locus of points from which the two tangents of the conic are
    perpendicular.

    Return value by conic type:

    - *Ellipses*: a real circle.
    - *Hyperbolas*: a real circle, zero-radius circle or imaginary circle for
      acute, rectangular and obtuse hyperbolas, respectively.
    - *Imaginary ellipses*: an imaginary circle.
    - *Parabolas*: the line pair of the
      [directrix](#parabola.parabola_directrix) and the ideal line.
    - *Finite point conics and line pairs crossing at a finite point*: a
      zero-radius circle at that point.
    - *Ideal point conics and line pairs crossing at an ideal point*: the
      double ideal line.
    - *Double lines and the zero matrix*: the zero matrix.

    *Definition*: <https://en.wikipedia.org/wiki/Director_circle><br>
    *Formula*:
    [research/construction/director_circle.py](../src/research/construction/director_circle.py)
    """
    a, _, c, d, e, f = conic_coefficients(conic.adjugate())
    return Matrix(
        [
            [-f, 0, d],
            [0, -f, e],
            [d, e, -(a + c)],
        ],
    )


#: The circle at the origin with radius 1.
UNIT_CIRCLE: ImmutableMatrix = ImmutableMatrix(circle(ORIGIN, r=1))

#: The circle at the origin with radius 𝑖.
IMAGINARY_UNIT_CIRCLE: ImmutableMatrix = ImmutableMatrix.eye(3)
