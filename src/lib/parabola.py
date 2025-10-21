from sympy import Expr, Matrix, sqrt

from lib.conic_classification import IsParabola
from lib.point import PointToXY


def _ParabolaDirectrixFromAdjugate(parabola_adjugate: Matrix) -> Matrix:
    """Computes the directrix of a parabola represented as the adjugate of a
    conic matrix. See [ParabolaDirectrix](#parabola.ParabolaDirectrix) for the
    details.
    """
    a, _, _, _, c, _, d, e, f = parabola_adjugate
    if f.is_zero is False:
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return Matrix([d, e, (a + c) / -2])


def ParabolaDirectrix(parabola: Matrix) -> Matrix:
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
    [research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)
    """
    return _ParabolaDirectrixFromAdjugate(parabola.adjugate())


def _ParabolaFocusFromAdjugate(parabola_adjugate: Matrix) -> Matrix:
    """Computes the focus of a parabola represented as the adjugate of a
    conic matrix. See [ParabolaFocus](#parabola.ParabolaFocus) for the details.
    """
    a, _, _, _, c, _, d, e, f = parabola_adjugate
    if f.is_zero is False:
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return parabola_adjugate * Matrix([d, e, -(a + c) / 2])


def ParabolaFocus(parabola: Matrix) -> Matrix:
    """Computes the focus of a parabola represented as a conic matrix.

    Special cases for other conic types:
    - Returns `[0, 0, 0]ᵀ` for
      - degenerate conics with one ideal point;
      - conics containing the ideal line.
    - Raises `ValueError` if the conic provably has 0 or 2 ideal points.
    - Returns an unspecified 3D column vector in all other cases.

    *Formula*:
    [research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)
    """
    return _ParabolaFocusFromAdjugate(parabola.adjugate())


def ParabolaVertex(parabola: Matrix) -> Matrix:
    """Computes the parabola's vertex.

    Returns a 2d vector.

    *Formula*: [research/parabola_vertex.py](../src/research/parabola_vertex.py)
    """
    a, _, _, b, c, _, d, e, _ = parabola
    focus_x, focus_y = PointToXY(ParabolaFocus(parabola))
    vertex_x = focus_x - (b * e - c * d) / (2 * (a + c) ** 2)
    vertex_y = focus_y - (b * d - a * e) / (2 * (a + c) ** 2)
    return Matrix([vertex_x, vertex_y])


def ParabolaDirection(parabola: Matrix) -> Matrix:
    """Computes the direction of a parabola modulo 2π.

    Unlike [FocalAxisDirection](#conic.FocalAxisDirection), which determines
    the direction only modulo π, this function resolves the full orientation.

    Returns the ideal point representing the parabola’s direction. This point
    coincides with both the ideal point lying on the parabola and its
    [ProjectiveConicCenter](#conic.ProjectiveConicCenter).
    """
    if IsParabola(parabola) is False:
        raise ValueError("Not a parabola")
    x, y, _ = parabola.row(0).cross(parabola.row(1))
    return Matrix([x, y, 0])


def ParabolaAxis(parabola: Matrix) -> Matrix:
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


def ParabolaFocalParameter(parabola: Matrix) -> Expr:
    """Computes the parabola's focus-directrix distance.

    *Formula*: [research/focal_parameter.py](../src/research/focal_parameter.py)
    """
    a, c, _ = parabola.diagonal()
    return sqrt(-parabola.det() / (a + c) ** 3)
