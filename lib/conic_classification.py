from sympy import Eq, Matrix
from sympy.core.logic import fuzzy_and


def IsCircular(conic: Matrix):
    """Whether there is a single center point around which the conic is
    invariant under all rotations.

    Circles, complex circles, zero-radius circles have such circular symmetry.
    Double ideal lines are not considered circular.
    """
    return fuzzy_and([conic[0].is_nonzero, conic[3].is_zero, Eq(conic[0], conic[4])])
