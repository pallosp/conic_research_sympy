from sympy import Matrix, pi, sqrt

from lib.ellipse import Ellipse
from lib.matrix import IsNonZeroMultiple
from lib.circle import UNIT_CIRCLE, Circle, CircleRadius, DirectorCircle


def test_unit_circle():
    assert IsNonZeroMultiple(UNIT_CIRCLE, Matrix.diag([1, 1, -1]))


def test_circle():
    circle = Circle((1, 2), 3)
    assert IsNonZeroMultiple(circle, Matrix([[-1, 0, 1], [0, -1, 2], [1, 2, 4]]))


def test_circle_radius():
    assert CircleRadius(Circle((1, 2), 3)) == 3
    assert CircleRadius(UNIT_CIRCLE * -2) == 1


def test_director_circle():
    assert DirectorCircle(Circle((1, 2), 3)) == Circle((1, 2), 3 * sqrt(2))
    assert DirectorCircle(Ellipse((1, 2), 3, 4)) == Circle((1, 2), 5)
    assert DirectorCircle(Ellipse((1, 2), 3, 4, r1_angle=pi / 4)) == Circle((1, 2), 5)
