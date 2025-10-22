from collections.abc import Sequence

from sympy import Matrix

from lib.matrix import is_nonzero_multiple


def are_projective_sets_equal(
    set1: Sequence[Matrix],
    set2: Sequence[Matrix],
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
