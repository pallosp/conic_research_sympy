from sympy import expand
from sympy import Matrix


def IsScalarMultiple(m1: Matrix, m2: Matrix) -> bool:
    """Tells whether two matrices are scalar multiples of each other."""
    if m1.shape != m2.shape:
        return False
    v1 = Matrix(list(m1))
    v2 = Matrix(list(m2))
    return expand(v1.dot(v1) * v2.dot(v2) - v1.dot(v2) ** 2) == 0
