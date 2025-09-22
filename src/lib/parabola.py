from sympy import Matrix


def ParabolaDirectrix(parabola: Matrix) -> Matrix:
    """Computes the directrix of a parabola represented as a conic matrix.

    Special cases for other conic types:
    - Returns the zero vector for coincident line pairs.
    - Returns the ideal line for non-coincident parallel line pairs.
    - Returns the ideal line for conics consisting of one finite and one ideal line.
    - Raises `ValueError` if the conic provably has 0 or 2 ideal points.
    - Returns an unspecified 3D column vector in all other cases.
    """
    a11, _, _, _, a22, _, a31, a32, a33 = parabola.adjugate()
    if a33.is_zero not in (True, None):
        raise ValueError("Not a parabola (or a degenerate conic with one ideal point)")
    return Matrix([a31, a32, (a11 + a22) / -2])
