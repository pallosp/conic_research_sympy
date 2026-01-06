from collections.abc import Callable

from sympy import Expr, Matrix, simplify
from sympy.core.logic import fuzzy_and


def is_homography(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Checks whether a transformation matrix represents a homography.

    A homography (projective transformation) matrix is a non-singular 3×3
    matrix, up to scale.

    The optional `simplifier` is applied to the determinant before comparing it
    to zero. If the determinant cannot be decided to be zero or non-zero after
    simplification, the function returns `None`.
    """
    if transformation.shape != (3, 3):
        return False

    return simplifier(transformation.det()).is_nonzero


def is_affine_transform(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Checks whether a matrix represents an affine transformation.

    An affine transformation is represented (up to scale) by a non-singular 3×3
    matrix whose last row is proportional to (0, 0, 1). Geometrically, affine
    transformations preserve parallelism and ratios of distances along a line,
    but not necessarily angles or lengths.

    The optional `simplifier` is applied to the affinity checking polynomials
    before comparing them to zero. If the polynomials cannot be decided to be
    zero or non-zero after simplification, the function returns `None`.
    """
    if transformation.shape != (3, 3):
        return False

    row2 = transformation.row(2)
    return fuzzy_and(
        [
            is_homography(transformation, simplifier=simplifier),
            simplifier(row2[0]).is_zero,
            simplifier(row2[1]).is_zero,
            simplifier(row2[2]).is_nonzero,
        ],
    )


def is_similarity(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Checks whether a matrix represents a similarity transformation.

    The optional `simplifier` is applied to the similarity checking polynomials
    before comparing them to zero. If the polynomials cannot be decided to be
    zero or non-zero after simplification, the function returns `None`.
    """
    if transformation.shape != (3, 3):
        return False

    a, b, c, d = transformation[:2, :2]

    return fuzzy_and(
        [
            is_affine_transform(transformation, simplifier=simplifier),
            simplifier(a * a - d * d).is_zero,
            simplifier(b * b - c * c).is_zero,
            simplifier(a * b + c * d).is_zero,
        ],
    )


def is_congruence(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Checks whether a matrix represents a congruence transformation.

    The optional `simplifier` is applied to the congruence checking polynomials
    before comparing them to zero. If the polynomials cannot be decided to be
    zero or non-zero after simplification, the function returns `None`.
    """
    if transformation.shape != (3, 3):
        return False

    a, b, *_, i = transformation

    return fuzzy_and(
        [
            is_similarity(transformation, simplifier=simplifier),
            simplifier(a * a + b * b - i * i).is_zero,
        ],
    )


def is_involution(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Checks whether applying the transformation twice is the identity.

    The optional `simplifier` is applied to the matrix elements before
    comparing them to zero. If they cannot be decided to be zero or non-zero
    after simplification, the function returns `None`.
    """
    if transformation.shape != (3, 3):
        return False

    double_transformation = simplifier(transformation * transformation)
    double_transformation /= double_transformation[0, 0]
    return (double_transformation - Matrix.eye(3)).is_zero_matrix
