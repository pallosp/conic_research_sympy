from collections.abc import Callable

from sympy import Expr, Matrix, simplify
from sympy.core.logic import fuzzy_and


def is_homography(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Tells whether a transformation matrix is a homography.

    Takes an optional `simplifier` callback that simplifies the determinant
    before it gets compared to zero. Returns `None` if the result is
    undecidable.
    """
    return fuzzy_and(
        [
            transformation.shape == (3, 3),
            simplifier(transformation.det()).is_nonzero,
        ],
    )


def is_affine_transform(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Tells whether a transformation matrix is an affine transformation.

    Takes an optional `simplifier` callback that simplifies the affinity
    checking polynomials before they get compared to zero. Returns `None` if
    undecidable.
    """
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
    """Tells whether a transformation matrix is a similarity transformation.

    Takes an optional `simplifier` callback that simplifies the similarity
    checking polynomials before they get compared to zero. Returns `None` if the
    result is undecidable.
    """
    if transformation.shape != (3, 3):
        return False

    a, b, _, c, d, _, _, _, _ = transformation

    return fuzzy_and(
        [
            is_affine_transform(transformation, simplifier=simplifier),
            simplifier(a**2 + c**2 - b**2 - d**2).is_zero,
            simplifier(a * b + c * d).is_zero,
        ],
    )
