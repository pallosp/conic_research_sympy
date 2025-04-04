from sympy import Matrix

from lib.matrix import IsScalarMultiple
from lib.circle import Circle


def test_circle():
    circle = Circle(0, 0, 1)
    assert IsScalarMultiple(circle, Matrix.diag([1, 1, -1]))
