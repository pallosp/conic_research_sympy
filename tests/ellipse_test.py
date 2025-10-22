import pytest
from sympy import Matrix, factor, nan, pi, simplify, symbols

from lib.central_conic import conic_center, semi_axis_lengths
from lib.circle import UNIT_CIRCLE, circle
from lib.conic_classification import is_ellipse
from lib.degenerate_conic import double_line_conic
from lib.ellipse import (
    ellipse,
    ellipse_from_foci_and_point,
    steiner_ellipse,
    steiner_inellipse,
)
from lib.incidence import conic_contains_point
from lib.line import line_between
from lib.matrix import is_nonzero_multiple
from lib.point import ORIGIN, centroid, point_to_vec3


class TestEllipseFromParams:
    def test_unit_circle_ellipse(self):
        assert ellipse(ORIGIN, 1, 1) == UNIT_CIRCLE

    def test_overspecified(self):
        with pytest.raises(ValueError, match="not both"):
            ellipse(ORIGIN, 1, 1, r1_angle=0, r1_direction=(1, 0))

    def test_wrong_r1_direction(self):
        with pytest.raises(ValueError, match="must be a 2D vector or an ideal point"):
            ellipse(ORIGIN, 1, 1, r1_direction=(1, 1, 1))

    def test_center(self):
        center = Matrix(symbols("x,y"))
        r1, r2 = symbols("r1,r2")
        r1_direction = symbols("r1_x,r1_y")
        symbolic_ellipse = ellipse(center, r1, r2, r1_direction=r1_direction)
        assert center == conic_center(symbolic_ellipse).applyfunc(factor)

    def test_axis_direction_vector_length_invariance(self):
        center = symbols("x,y")
        r1, r2, r1_dir_x, r1_dir_y = symbols("r1,r2,dx,dy", positive=True)
        conic1 = ellipse(center, r1, r2, r1_direction=(r1_dir_x, r1_dir_y))
        conic2 = ellipse(center, r1, r2, r1_direction=(r1_dir_x * -2, r1_dir_y * -2))
        assert is_nonzero_multiple(conic1, conic2)

    def test_axis_lengths_circle(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = ellipse(center, r, r)
        assert semi_axis_lengths(circle) == (r, r)

    def test_axis_lengths_numeric(self):
        ellipse1 = ellipse((1, 2), 3, 4, r1_direction=(5, 6))
        assert semi_axis_lengths(ellipse1) == (3, 4)
        ellipse2 = ellipse((1, 2), 3, 4, r1_angle=pi / 6)
        assert semi_axis_lengths(ellipse2) == (3, 4)

    def test_axis_lengths_general_case(self):
        center = symbols("x,y")
        r_min, r_diff = symbols("r_min,r_diff", positive=True)
        r_max = r_min + r_diff
        symbolic_ellipse = ellipse(center, r_min, r_max, r1_direction=(73, -25))
        axes = [factor(simplify(leng)) for leng in semi_axis_lengths(symbolic_ellipse)]
        assert axes in ([r_min, r_max], [r_max, r_min])


class TestEllipseFromFociAndPoint:
    def test_numeric_ellipse(self):
        f1 = (1, 2)
        f2 = (3, 4)
        p = (0, 0)
        ellipse = ellipse_from_foci_and_point(f1, f2, p)
        assert is_ellipse(ellipse)
        center = [coord.simplify() for coord in conic_center(ellipse)]
        assert center == [2, 3]
        assert conic_contains_point(ellipse, p)

    def test_collinear_points(self):
        f1 = (1, 2)
        f2 = (3, 4)
        conic = ellipse_from_foci_and_point(f1, f2, (2, 3))
        assert is_nonzero_multiple(conic, double_line_conic(line_between(f1, f2)))
        conic = ellipse_from_foci_and_point(f1, f2, (1, 2))
        assert is_nonzero_multiple(conic, double_line_conic(line_between(f1, f2)))

    def test_coincident_foci(self):
        f = (2, 3)
        p = (5, 7)
        circular_ellipse = ellipse_from_foci_and_point(f, f, p)
        assert is_nonzero_multiple(circular_ellipse, circle(f, 5))

    def test_coincident_foci_and_point(self):
        assert ellipse_from_foci_and_point((1, 2), (1, 2), (1, 2)).is_zero_matrix

    def test_ideal_point_focus(self):
        f1 = (1, 2)
        f2 = (3, 4, 0)
        p = (5, 6)
        assert nan in ellipse_from_foci_and_point(f1, f2, p)

    def test_ideal_incident_point(self):
        f1 = (1, 2)
        f2 = (3, 4)
        p = (5, 6, 0)
        assert nan in ellipse_from_foci_and_point(f1, f2, p)


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
