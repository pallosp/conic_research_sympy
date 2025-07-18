from sympy import Eq, Matrix
from sympy.core.logic import fuzzy_and


def IsDegenerate(conic: Matrix) -> bool | None:
    return conic.det().is_zero


def IsNonDegenerate(conic: Matrix) -> bool | None:
    return conic.det().is_nonzero


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
