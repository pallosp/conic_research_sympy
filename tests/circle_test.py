from sympy import Matrix, pi, sqrt
from sympy.abc import x, y

from lib.circle import UNIT_CIRCLE, circle, circle_radius, director_circle
from lib.conic import conic_from_poly
from lib.ellipse import ellipse
from lib.matrix import is_nonzero_multiple


def test_unit_circle():
    assert is_nonzero_multiple(UNIT_CIRCLE, Matrix.diag([1, 1, -1]))


def test_circle():
    circle_matrix = circle((1, 2), 3)
    expected = Matrix([[-1, 0, 1], [0, -1, 2], [1, 2, 4]])
    assert is_nonzero_multiple(circle_matrix, expected)


def test_circle_radius():
    assert circle_radius(circle((1, 2), 3)) == 3
    assert circle_radius(UNIT_CIRCLE * -2) == 1


class TestDirectorCircle:
    def test_axis_aligned_ellipse(self):
        assert director_circle(circle((1, 2), 3)) == circle((1, 2), 3 * sqrt(2))
        assert director_circle(ellipse((1, 2), 3, 4)) == circle((1, 2), 5)

    def test_rotated_ellipse(self):
        rotated_ellipse = ellipse((1, 2), 3, 4, r1_angle=pi / 4)
        assert director_circle(rotated_ellipse) == circle((1, 2), 5)

    def test_rectangular_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 1)
        assert director_circle(hyperbola) == circle((0, 0), 0)
