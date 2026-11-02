from sympy import ImmutableMatrix, MatrixExpr

from lib.matrix import (
    outer_product_root,
    pivot_column,
    pivot_column_and_row,
    skew_matrix,
)
from lib.point import point_to_vec3
from lib.type_aliases import Conic, Line, Point, PointVector

#: The zero matrix. Every point satisfies its equation.
ZERO_CONIC: ImmutableMatrix = ImmutableMatrix.zeros(3, 3)


def line_pair_conic(line1: Line, line2: Line) -> Conic:
    """Constructs a degenerate conic consisting of two lines.

    Formula: `C = (l₁ l₂ᵀ + l₂ l₁ᵀ) / 2`

    Derivation: a point `p` lies on line `l` iff `pᵀ l = lᵀ p = 0`, so it
    lies on the line pair iff

        (pᵀ l₁)(l₂ᵀ p) = pᵀ (l₁ l₂ᵀ) p = 0

    The matrix `M = l₁ l₂ᵀ` defines this quadratic form but is not symmetric.
    Since `pᵀ M p` is a number, it equals its transpose `pᵀ Mᵀ p`. Averaging
    the two gives the same quadratic form with the symmetric matrix
    `(M + Mᵀ) / 2`.
    """
    if line1.shape != (3, 1) or line2.shape != (3, 1):
        raise ValueError("The lines must be 3-dimensional column vectors.")
    return (line1 * line2.T + line2 * line1.T) / 2


def double_line_conic(line: Line) -> Conic:
    """Constructs a degenerate conic consisting of two coincident lines."""
    if line.shape != (3, 1):
        raise ValueError("The line must be a 3-dimensional column vector.")
    return line * line.T


def point_conic(point: Point) -> Conic:
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

    An alternative way to construct the point conic is to take the tangent lines
    through the point to an imaginary ellipse. The choice of the diag(-1, -1, -1)
    imaginary unit circle yields the same conic matrix.

    Hint: Use [extract_point](#degenerate_conic.extract_point) to recover the
    point from the resulting conic.
    """
    skew = skew_matrix(point_to_vec3(point))
    return skew * skew  # = skew.T * Matrix.diag([-1, -1, -1]) * skew


def split_to_lines(
    conic: Conic,
) -> tuple[Line | MatrixExpr, Line | MatrixExpr]:
    """Splits a degenerate conic into two lines.

    Special cases:
     - For non-degenerate conics the result is unspecified.
     - For point conics the lines are complex conjugates.
     - Returns two `[0, 0, 0]ᵀ` vectors for the zero matrix.
     - For symbolic conics the lines may be unevaluated `MatrixExpr`s, or
       explicit matrices with unevaluated elements such as
       `OuterProductRoot(...)[0, 0]`. Substitutions evaluate them once the
       result becomes decidable.

    *Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 11.1

    The conic is `C = l₁·l₂ᵀ + l₂·l₁ᵀ` for some lines `l₁` and `l₂` (see
    [line_pair_conic](#degenerate_conic.line_pair_conic)). Let `q = l₁⨯l₂`
    be their intersection point.

     1. `C - skew(q) = 2·l₁·l₂ᵀ` and `C + skew(q) = 2·l₂·l₁ᵀ`, because
        `skew(l₁⨯l₂) = l₂·l₁ᵀ - l₁·l₂ᵀ` (see [skew_matrix](#matrix.skew_matrix)).
     2. `adj(C) = -q·qᵀ`
        ([proof](../src/research/matrices/line_pair_adjugate.py)), so
        [outer_product_root](#matrix.outer_product_root) recovers `q` from
        `C` up to sign. The sign only swaps the lines. For point conics `q`
        is imaginary; for double lines it's zero.
     3. `l₁·l₂ᵀ` has rank 1: its columns are multiples of `l₁` and its rows
        are multiples of `l₂`.
        [pivot_column_and_row](#matrix.pivot_column_and_row) takes the column
        and the row of a non-zero element.
    """
    q = outer_product_root(-conic.adjugate())
    # SymPy's stubs reject Matrix + MatrixExpr (fix upstream).
    return pivot_column_and_row(conic + skew_matrix(q))  # pyright: ignore[reportOperatorIssue]


def extract_point(degenerate_conic: Conic) -> PointVector | MatrixExpr:
    """Extracts the point from a point conic or the intersection of the
    lines from a line pair conic.

    This point is the conic's *singular point*: the point `p` for which
    `C·p = 0`, i.e. whose [polar line](#conic.polar_line) is undefined.
    Non-degenerate conics have no singular points, while every point of a
    double line is singular.

    Special cases:
     - Returns `[0, 0, 0]ᵀ` for double line pairs.
     - Returns an unspecified 3D column vector if the conic is not degenerate.
     - May return an unevaluated `MatrixExpr` for symbolic conics.

    *Algorithm*:

    Let the rows of the conic matrix be `r₁`, `r₂` and `r₃`. The point `p` is
    on the conic iff `pᵀ C p = 0`, i.e. `p·(r₁·p, r₂·p, r₃·p) = 0`.

    `p = r₁⨯r₂` is a solution, because `r₁·(r₁⨯r₂) = 0`, `r₂·(r₁⨯r₂) = 0` and
    `r₃·(r₁⨯r₂) = det C = 0`. So are `p = r₂⨯r₃` and `p = r₃⨯r₁`. When the
    conic matrix is a rank 2 matrix (point conic or non-coincident line pair),
    at least one of these is a non-zero vector.

    Since `C` is symmetric, these cross products are the columns of `adj(C)`.
    """
    return pivot_column(degenerate_conic.adjugate())
