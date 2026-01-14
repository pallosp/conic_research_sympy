from sympy import Expr, Matrix, sqrt

from lib.conic_classes import is_parabola
from lib.point import point_to_xy


def _parabola_directrix_from_adjugate(parabola_adjugate: Matrix) -> Matrix:
    """Computes the directrix of a parabola represented as the adjugate of a
    conic matrix. See [parabola_directrix](#parabola.parabola_directrix) for the
    details.
    """
    a, _, _, _, c, _, d, e, f = parabola_adjugate
    if f.is_zero is False:
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return Matrix([d, e, (a + c) / -2])


def parabola_directrix(parabola: Matrix) -> Matrix:
    """Computes the directrix of a parabola represented as a conic matrix.

    Special cases for other conic types:
    - Returns `[0, 0, 0]ᵀ` for coincident line pairs.
    - Returns the ideal line for
      - non-coincident parallel line pairs;
      - conics consisting of one finite and one ideal line;
      - ideal point conics.
    - Raises `ValueError` if the conic provably has 0 or 2 ideal points.
    - Returns an unspecified 3D column vector in all other cases.

    *Formula*:
    [research/conic_properties/focus_directrix_eccentricity.py](../src/research/conic_properties/focus_directrix_eccentricity.py)
    """
    return _parabola_directrix_from_adjugate(parabola.adjugate())


def _parabola_focus_from_adjugate(parabola_adjugate: Matrix) -> Matrix:
    """Computes the focus of a parabola represented as the adjugate of a
    conic matrix. See [parabola_focus](#parabola.parabola_focus) for the details.
    """
    a, _, _, _, c, _, d, e, f = parabola_adjugate
    if f.is_zero is False:
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return parabola_adjugate * Matrix([d, e, -(a + c) / 2])


def parabola_focus(parabola: Matrix) -> Matrix:
    """Computes the focus of a parabola represented as a conic matrix.

    Special cases for other conic types:
    - Returns `[0, 0, 0]ᵀ` for
      - degenerate conics with one ideal point;
      - conics containing the ideal line.
    - Raises `ValueError` if the conic provably has 0 or 2 ideal points.
    - Returns an unspecified 3D column vector in all other cases.

    *Formula*:
    [research/conic_properties/focus_directrix_eccentricity.py](../src/research/conic_properties/focus_directrix_eccentricity.py)
    """
    return _parabola_focus_from_adjugate(parabola.adjugate())


def parabola_vertex(parabola: Matrix) -> Matrix:
    """Computes the parabola's vertex.

    Returns the point's coordinates as a 2D column vector.

    *Formula*:
    [research/conic_properties/parabola_vertex.py](../src/research/conic_properties/parabola_vertex.py)
    """
    a, _, _, b, c, _, d, e, _ = parabola
    focus_x, focus_y = point_to_xy(parabola_focus(parabola))
    vertex_x = focus_x - (b * e - c * d) / (2 * (a + c) ** 2)
    vertex_y = focus_y - (b * d - a * e) / (2 * (a + c) ** 2)
    return Matrix([vertex_x, vertex_y])


def parabola_direction(parabola: Matrix) -> Matrix:
    """Computes the direction of a parabola modulo 2π.

    Unlike [focal_axis_direction](#conic_direction.focal_axis_direction), which
    determines the direction only modulo π, this function resolves the full
    orientation.

    Returns the ideal point representing the parabola’s direction. This point
    coincides with both the ideal point lying on the parabola and its
    [projective_conic_center](#conic.projective_conic_center).
    """
    if is_parabola(parabola) is False:
        raise ValueError("Not a parabola")
    x, y, _ = parabola.row(0).cross(parabola.row(1))
    return Matrix([x, y, 0])


def parabola_axis(parabola: Matrix) -> Matrix:
    """Computes the parabola's focal axis line.

    It's the polar line corresponding to the ideal point on the directrix.

    Special cases for other conic types:
    - Returns `[0, 0, 0]ᵀ` for
      - degenerate conics with one ideal point;
      - conics containing the ideal line.
    - Raises `ValueError` if the conic provably has 0 or 2 ideal points.
    - Returns an unspecified 3D column vector in all other cases.
    """
    x, y, discriminant = parabola.row(0).cross(parabola.row(1))
    if discriminant.is_zero is False:
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return parabola * Matrix([-y, x, 0])


def parabola_focal_parameter(parabola: Matrix) -> Expr:
    """Computes the parabola's focus-directrix distance.

    *Formula*:
    [research/conic_properties/focal_parameter.py](../src/research/conic_properties/focal_parameter.py)
    """
    a, c, _ = parabola.diagonal()
    return sqrt(-parabola.det() / (a + c) ** 3)
