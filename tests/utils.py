from collections.abc import Sequence
from typing import TypeGuard

from sympy import MatrixBase

from lib.matrix import is_nonzero_multiple


def are_projective_sets_equal(
    set1: Sequence[MatrixBase],
    set2: Sequence[MatrixBase],
) -> bool:
    """Compares two sets of projective points or lines for equality."""
    if len(set1) != len(set2):
        return False

    remaining = list(set2)

    for v1 in set1:
        found = False
        for v2 in remaining:
            if is_nonzero_multiple(v1, v2):
                remaining.remove(v2)
                found = True
                break
        if not found:
            return False
    return True


def is_explicit_matrix(matrix: object) -> TypeGuard[MatrixBase]:
    """Tells whether a matrix is explicit, as opposed to an unevaluated
    `MatrixExpr` such as a `MatrixFunction`.
    """
    return isinstance(matrix, MatrixBase)
