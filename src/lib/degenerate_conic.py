from collections.abc import Sequence

from sympy import Expr, Function, Matrix, sqrt

from lib.matrix import ConicMatrix, NonZeroCross, SkewMatrix
from lib.point import PointToVec3


def LinePair(line1: Matrix, line2: Matrix) -> Matrix:
    """Constructs a conic section from two projective lines."""
    if line1.shape != (3, 1) or line2.shape != (3, 1):
        raise ValueError("The lines must be 3-dimensional column vectors.")
    return (line1 * line2.T + line2 * line1.T) / 2


def PointConic(point: Matrix | Sequence[Expr]) -> Matrix:
    """Constructs a conic that degenerates to a single point.

    Let the point's homogeneous coordinates be `(x, y, z)` and let
    the variable point on the conic be `v = (X, Y, Z)`. The conic is
    defined by the quadratic form `vᵀ C v = 0`.

    For a conic consisting of just the given point, the condition
    `x : y : z = X : Y : Z` must hold. Such conics can be expressed by the
    equation

        λ₁(Y*z - Z*y)² + λ₂(Z*x - X*z)² + λ₃(X*y - Y*x)² = 0

    where `λ₁`, `λ₂` and `λ₃` are either all positive or all negative.
    This function returns the corresponding conic matrix for the choice
    `λ₁ = λ₂ = λ₃ = -1`.
    """
    x, y, z = PointToVec3(point)
    return ConicMatrix(
        -y * y - z * z,
        x * y,
        -z * z - x * x,
        x * z,
        y * z,
        -x * x - y * y,
    )


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
        a, c, f = adj.diagonal()

        # Lemma: If a symmetric 3x3 matrix is singular, the diagonal elements
        # of its adjugate matrix (a, c, f) are either all ≥0 or all ≤0.
        if a.is_nonzero:
            conic = conic + SkewMatrix(adj.col(0) / sqrt(-a))
        elif c.is_nonzero:
            conic = conic + SkewMatrix(adj.col(1) / sqrt(-c))
        elif f.is_nonzero:
            conic = conic + SkewMatrix(adj.col(2) / sqrt(-f))
        elif not (a.is_zero and c.is_zero and f.is_zero):
            return None

        cross = NonZeroCross(conic)
        if isinstance(cross, NonZeroCross):
            return None
        return (cross[0], cross[1].T)
