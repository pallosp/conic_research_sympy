from sympy import Matrix
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
