from collections.abc import Callable

from sympy import Expr, Matrix, simplify
from sympy.core.logic import fuzzy_and


def is_homography(transformation: Matrix) -> bool | None:
    """Tells whether a transformation matrix is a homography.

    Returns None if undecidable.
    """
    return fuzzy_and(
        [
            transformation.shape == (3, 3),
            transformation.det().simplify().is_nonzero,
        ],
    )


def is_affine_transform(transformation: Matrix) -> bool | None:
    """Tells whether a transformation matrix is an affine transformation.

    Returns None if undecidable.
    """
    row2 = transformation.row(2)
    return fuzzy_and(
        [
            is_homography(transformation),
            row2[0].is_zero,
            row2[1].is_zero,
            row2[2].is_nonzero,
        ],
    )


def is_similarity(
    transformation: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = simplify,
) -> bool | None:
    """Tells whether a transformation matrix is a similarity transformation.

    Takes an optional `simplifier` callback that simplifies the similarity
    polynomials before they get compared to zero. Returns `None` if the result
    is undecidable.
    """
    if transformation.shape != (3, 3):
        return False

    a, b, _, c, d, _, _, _, _ = transformation

    return fuzzy_and(
        [
            is_affine_transform(transformation),
            simplifier(a**2 + c**2 - b**2 - d**2).is_zero,
            simplifier(a * b + c * d).is_zero,
        ],
    )
