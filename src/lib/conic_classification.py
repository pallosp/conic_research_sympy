from sympy import Eq, Matrix
from sympy.core.logic import fuzzy_and, fuzzy_not

from lib.matrix import IsDefinite


def IsDegenerate(conic: Matrix) -> bool | None:
    """Tells whether the conic is degenerate.

    Degenerate conics consist of a single projective point or a pair of
    projective lines. The zero matrix is also considered degenerate.
    Returns None if undecidable.
    """
    return conic.det().is_zero


def IsNonDegenerate(conic: Matrix) -> bool | None:
    """Tells whether the conic is non-degenerate.

    Non-degenerate conics include real or complex ellipses, parabolas and
    hyperbolas. Returns None if undecidable.
    """
    return conic.det().is_nonzero


def IsFiniteConic(conic: Matrix) -> bool | None:
    """Tells whether all points on the conic are finite.

    Returns None if undecidable.
    """
    return conic[:2, :2].det().factor().is_positive


def IsComplexEllipse(conic: Matrix) -> bool | None:
    """Tells Whether the conic is an ellipse with a real center and imaginary
    radii.

    Returns None if undecidable.
    """
    return IsDefinite(conic)


def IsEllipse(conic: Matrix) -> bool | None:
    """Tells whether the conic is an ellipse with real radii.

    Returns None if undecidable.
    """
    return fuzzy_and(
        [
            IsNonDegenerate(conic),
            conic[:2, :2].det().is_positive,
            fuzzy_not(IsDefinite(conic)),
        ]
    )


def IsParabola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a parabola.

    Returns None if undecidable.
    """
    return fuzzy_and([IsNonDegenerate(conic), conic[:2, :2].det().is_zero])


def IsHyperbola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a hyperbola.

    Returns None if undecidable.
    """
    return fuzzy_and([IsNonDegenerate(conic), conic[:2, :2].det().is_negative])


def IsCircular(conic: Matrix) -> bool | None:
    """Tells whether there is a single center point around which the conic is
    invariant under all rotations.

    Circles, complex circles, zero-radius circles have such circular symmetry.
    Double ideal lines are not considered circular. Returns None if undecidable.
    """
    return fuzzy_and([conic[0].is_nonzero, conic[3].is_zero, Eq(conic[0], conic[4])])


def IsLinePair(conic: Matrix) -> bool | None:
    """Tells whether the conic is the union of two projective lines.

    Returns None if undecidable.
    """
    return fuzzy_and(
        [
            IsDegenerate(conic),
            fuzzy_not(IsFiniteConic(conic)),
            fuzzy_not(conic.is_zero_matrix),
        ]
    )


def IsDoubleLine(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of two coincident projective lines.

    Returns None if undecidable.
    """
    # The conic represents a double line iff its matrix has rank 1.
    #
    # A matrix has rank ≤ 1 iff the cross product of every pair of columns is
    # zero, or equivalently its adjugate is the zero matrix.
    return fuzzy_and(
        [
            conic.adjugate().is_zero_matrix,
            fuzzy_not(conic.is_zero_matrix),
        ]
    )
