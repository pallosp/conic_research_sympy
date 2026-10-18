from sympy import ImmutableMatrix, Matrix, MatrixBase, MatrixExpr

from lib.matrix import (
    conic_matrix,
    outer_product_root,
    pivot_column,
    pivot_column_and_row,
    skew_matrix,
)
from lib.point import point_to_vec3
from lib.type_aliases import Point

#: The zero matrix. Every point satisfies its equation.
ZERO_CONIC: ImmutableMatrix = ImmutableMatrix.zeros(3, 3)


def line_pair_conic(line1: MatrixBase, line2: MatrixBase) -> Matrix:
    """Constructs a conic section from two projective lines."""
    if line1.shape != (3, 1) or line2.shape != (3, 1):
        raise ValueError("The lines must be 3-dimensional column vectors.")
    return (line1 * line2.T + line2 * line1.T) / 2


def double_line_conic(line: MatrixBase) -> MatrixBase:
    """Constructs a degenerate conic consisting of two coincident lines."""
    if line.shape != (3, 1):
        raise ValueError("The line must be a 3-dimensional column vector.")
    return line * line.T


def point_conic(point: Point) -> Matrix:
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

    Hint: Use [extract_point](#degenerate_conic.extract_point) to recover the
    point from the resulting conic.
    """
    x, y, z = point_to_vec3(point)
    return conic_matrix(
        -y * y - z * z,
        x * y,
        -z * z - x * x,
        x * z,
        y * z,
        -x * x - y * y,
    )


def split_to_lines(
    conic: MatrixBase,
) -> tuple[MatrixBase | MatrixExpr, MatrixBase | MatrixExpr]:
    """Splits a degenerate conic into two lines.

    Special cases:
     - For non-degenerate conics the result is unspecified.
     - For point conics the lines will be complex conjugates.
     - Returns two `[0, 0, 0]ᵀ` vectors for the zero matrix.
     - For symbolic conics the lines may be unevaluated `MatrixExpr`s, or
       explicit matrices with unevaluated elements such as
       `OuterProductRoot(...)[0, 0]`. Substitutions evaluate them once the
       result becomes decidable.

    *Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 11.1
    """
    # The conic consisting of the lines l and m is represented by the
    # symmetric matrix C = (l·mᵀ + m·lᵀ)/2, or a multiple of it.
    #
    # Goal: find the rank-1 matrix l·mᵀ. All its columns are multiples of l
    # and all its rows are multiples of m, so any non-zero element yields
    # both lines. pivot_column_and_row does this final step.
    #
    # l·mᵀ differs from C by an antisymmetric matrix, which can be written
    # as the skew matrix of a vector q:
    #
    #   l·mᵀ = C + (l·mᵀ - m·lᵀ)/2 = C + skew(q),   where q = (m⨯l)/2
    #
    # q is a multiple of the lines' intersection point. It's unknown, but
    # it can be computed from C alone, because adj(C) = -q·qᵀ.
    # outer_product_root recovers q from q·qᵀ up to sign. The sign doesn't
    # matter: C - skew(q) = m·lᵀ, which yields the same lines in swapped
    # order.
    #
    # Special cases:
    #  - Point conics: l and m are complex conjugates, so q is imaginary.
    #  - Double lines: l = m, so q = 0 and C = l·lᵀ is split directly.
    #  - Symbolic conics: outer_product_root and pivot_column_and_row return
    #    unevaluated matrix expressions when they can't decide which element
    #    is non-zero. Each step accepts such expressions, so the lines
    #    evaluate once substitutions make every step decidable.
    q = outer_product_root(-conic.adjugate())
    return pivot_column_and_row(conic + skew_matrix(q))


def extract_point(degenerate_conic: MatrixBase) -> MatrixBase | MatrixExpr:
    """Extracts the point from a point conic or the intersection of the
    lines from a line pair conic.

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
