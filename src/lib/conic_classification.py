from sympy import Eq, Matrix
from sympy.core.logic import fuzzy_and, fuzzy_not

from src.lib.matrix import IsDefinite


def IsDegenerate(conic: Matrix) -> bool | None:
    return conic.det().is_zero


def IsNonDegenerate(conic: Matrix) -> bool | None:
    return conic.det().is_nonzero


def IsFiniteConic(conic: Matrix) -> bool | None:
    """Whether all points on the conic are finite."""
    return conic[:2, :2].det().factor().is_positive


def IsComplexEllipse(conic: Matrix) -> bool | None:
    """Whether the conic is an ellipse with a real center and imaginary radii."""
    return IsDefinite(conic)


def IsEllipse(conic: Matrix) -> bool | None:
    """Whether the conic is an ellipse with real radii."""
    return fuzzy_and(
        [
            IsNonDegenerate(conic),
            conic[:2, :2].det().is_positive,
            fuzzy_not(IsDefinite(conic)),
        ]
    )


def IsParabola(conic: Matrix) -> bool | None:
    return fuzzy_and([IsNonDegenerate(conic), conic[:2, :2].det().is_zero])


def IsHyperbola(conic: Matrix) -> bool | None:
    return fuzzy_and([IsNonDegenerate(conic), conic[:2, :2].det() < 0])


def IsCircular(conic: Matrix) -> bool | None:
    """Whether there is a single center point around which the conic is
    invariant under all rotations.

    Circles, complex circles, zero-radius circles have such circular symmetry.
    Double ideal lines are not considered circular.
    """
    return fuzzy_and([conic[0].is_nonzero, conic[3].is_zero, Eq(conic[0], conic[4])])


def IsLinePair(conic: Matrix) -> bool | None:
    """Whether the conic is the union of two projective lines."""
    return fuzzy_and(
        [
            IsDegenerate(conic),
            fuzzy_not(IsFiniteConic(conic)),
            fuzzy_not(conic.is_zero_matrix),
        ]
    )
