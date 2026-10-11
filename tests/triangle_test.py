from sympy import Matrix, symbols
from sympy.abc import x, y

from lib.central_conic import conic_center
from lib.circle import circle
from lib.degenerate_conic import line_pair_conic
from lib.incidence import conic_contains_point
from lib.line import IDEAL_LINE, line_between
from lib.matrix import is_nonzero_multiple, is_positive_multiple
from lib.point import centroid, ideal_point, point_to_vec3
from lib.triangle import (
    circumcenter,
    circumcircle,
    steiner_ellipse,
    steiner_inellipse,
)


class TestCircumcircle:
    def test_homogeneous_points(self):
        assert is_nonzero_multiple(
            circumcircle((8, 4, 2), (-1, -5, -1), (-2, 2)), circle((1, 2), r=3)
        )

    def test_sign_depends_on_orientation(self):
        p1, p2, p3 = (4, 2), (1, 5), (-2, 2)  # counterclockwise
        expected = circle((1, 2), r=3)
        assert is_positive_multiple(circumcircle(p1, p2, p3), expected)
        assert is_positive_multiple(circumcircle(p2, p1, p3), -expected)

    def test_symbolic_points(self):
        p1, p2, p3 = (2, -1), (x, y), (6, 5)
        conic = circumcircle(p1, p2, p3)
        assert conic[0] == conic[4]
        assert conic[1] == 0
        assert conic_contains_point(conic, p1)
        assert conic_contains_point(conic, (x, y))
        assert conic_contains_point(conic, p3)

    def test_collinear_points(self):
        assert is_nonzero_multiple(
            circumcircle((0, 0), (1, 1), (3, 3)),
            line_pair_conic(line_between((0, 0), (1, 1)), IDEAL_LINE),
        )

    def test_one_ideal_point(self):
        assert is_nonzero_multiple(
            circumcircle((1, 2), (3, 4), (1, 0, 0)),
            line_pair_conic(line_between((1, 2), (3, 4)), IDEAL_LINE),
        )

    def test_two_ideal_points(self):
        assert circumcircle((1, 2), (1, 0, 0), (0, 1, 0)).is_zero_matrix

    def test_coincident_points(self):
        assert circumcircle((1, 2), (1, 2), (5, 7)).is_zero_matrix


class TestCircumcenter:
    def test_finite_points(self):
        assert is_nonzero_multiple(
            circumcenter((8, 4, 2), (-1, -5, -1), (-2, 2)), Matrix([1, 2, 1])
        )

    def test_symbolic_points(self):
        a, b, r = symbols("a b r", positive=True)
        points = [(a + r, b), (a, b + r), (a - r, b)]
        assert is_nonzero_multiple(circumcenter(*points), Matrix([a, b, 1]))

    def test_collinear_points(self):
        assert is_nonzero_multiple(
            circumcenter((0, 0), (1, 1), (3, 3)), ideal_point(1, -1)
        )

    def test_one_ideal_point(self):
        assert is_nonzero_multiple(
            circumcenter((1, 2), (3, 2), (1, 0, 0)), ideal_point(0, 1)
        )

    def test_two_ideal_points(self):
        assert circumcenter((1, 2), (1, 0, 0), (0, 1, 0)).is_zero_matrix

    def test_coincident_points(self):
        assert circumcenter((1, 2), (1, 2), (5, 7)).is_zero_matrix


class TestSteinerEllipse:
    def test_circumellipse(self):
        p1 = (4, 1)
        p2 = (7, 3)
        p3 = (5, 5)
        ellipse = steiner_ellipse(p1, p2, p3)
        assert conic_contains_point(ellipse, point_to_vec3(p1))
        assert conic_contains_point(ellipse, point_to_vec3(p2))
        assert conic_contains_point(ellipse, point_to_vec3(p3))
        assert conic_center(ellipse) == centroid(p1, p2, p3)

    def test_inellipse(self):
        p1 = (4, 1)
        p2 = (7, 3)
        p3 = (5, 5)
        ellipse = steiner_inellipse(p1, p2, p3)
        assert conic_contains_point(ellipse, point_to_vec3(centroid(p2, p3)))
        assert conic_contains_point(ellipse, point_to_vec3(centroid(p3, p1)))
        assert conic_contains_point(ellipse, point_to_vec3(centroid(p1, p2)))
        assert conic_center(ellipse) == centroid(p1, p2, p3)
