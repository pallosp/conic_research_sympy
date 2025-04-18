from sympy import Matrix

from lib.matrix import IsScalarMultiple
from lib.circle import Circle, UNIT_CIRCLE


def test_unit_circle():
    assert IsScalarMultiple(UNIT_CIRCLE, Matrix.diag([1, 1, -1]))


def test_circle():
    circle = Circle((1, 2), 3)
    assert IsScalarMultiple(circle, Matrix([[-1, 0, 1], [0, -1, 2], [1, 2, 4]]))
