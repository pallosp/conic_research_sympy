from sympy import Matrix, pi, sqrt
from sympy.abc import x, y

from lib.conic import ConicFromPoly
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


class TestDirectorCircle:
    def test_axis_aligned_ellipse(self):
        assert DirectorCircle(Circle((1, 2), 3)) == Circle((1, 2), 3 * sqrt(2))
        assert DirectorCircle(Ellipse((1, 2), 3, 4)) == Circle((1, 2), 5)

    def test_rotated_ellipse(self):
        ellipse = Ellipse((1, 2), 3, 4, r1_angle=pi / 4)
        assert DirectorCircle(ellipse) == Circle((1, 2), 5)

    def test_rectangular_hyperbola(self):
        hyperbola = ConicFromPoly(x * y - 1)
        assert DirectorCircle(hyperbola) == Circle((0, 0), 0)
