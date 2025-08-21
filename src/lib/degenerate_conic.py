from sympy import Function, Matrix, sqrt

from lib.matrix import NonZeroCross, SkewMatrix


def LinePair(line1: Matrix, line2: Matrix) -> Matrix:
    """Constructs a conic section from two projective lines."""
    assert line1.shape == (3, 1)
    assert line2.shape == (3, 1)
    return (line1 * line2.T + line2 * line1.T) / 2


class SplitToLines(Function):
    """Splits a degenerate conic into two lines.

    Special cases:
     - For non-degenerate conics the result is unspecified.
     - For point conics the lines will be complex conjugates.
     - For symbolic conics returns an unevaluated `sympy.Function`.

    Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 11.1
    """

    @classmethod
    def eval(cls, conic: Matrix) -> tuple[Matrix, Matrix] | None:
        adj = conic.adjugate()
        A, C, F = adj.diagonal()

        if A.is_nonzero:
            conic = conic + SkewMatrix(adj.col(0) / sqrt(-A))
        elif C.is_nonzero:
            conic = conic + SkewMatrix(adj.col(1) / sqrt(-C))
        elif F.is_nonzero:
            conic = conic + SkewMatrix(adj.col(2) / sqrt(-F))
        elif not (A.is_zero and C.is_zero and F.is_zero):
            return None

        cross = NonZeroCross(conic)
        if isinstance(cross, NonZeroCross):
            return None
        return (cross[0], cross[1].T)
