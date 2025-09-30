from sympy import Matrix


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
    a, _, _, _, c, _, d, e, f = parabola.adjugate()
    if f.is_zero not in (True, None):
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return Matrix([d, e, (a + c) / -2])


def ParabolaFocus(parabola: Matrix) -> Matrix:
    """Computes the focus of a parabola represented as a conic matrix.

    Special cases for other conic types:
    - Returns `[0, 0, 0]ᵀ` for
      - conics with one ideal point;
      - conics containing the ideal line.
    - Raises `ValueError` if the conic provably has 0 or 2 ideal points.
    - Returns an unspecified 3D column vector in all other cases.

    *Formula*:
    [research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)
    """
    adj = parabola.adjugate()
    a, _, _, _, c, _, d, e, f = adj
    if f.is_zero not in (True, None):
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return adj * Matrix([d, e, -(a + c) / 2])
