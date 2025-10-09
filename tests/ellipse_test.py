import pytest
from sympy import factor, nan, pi, simplify, symbols

from lib.central_conic import ConicCenter, SemiAxisLengths
from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import ConicContainsPoint
from lib.conic_classification import IsEllipse
from lib.degenerate_conic import DoubleLine
from lib.ellipse import (
    Ellipse,
    EllipseFromFociAndPoint,
    SteinerEllipse,
    SteinerInellipse,
)
from lib.line import LineBetween
from lib.matrix import IsNonZeroMultiple
from lib.point import ORIGIN, Centroid, PointToVec3


class TestEllipseFromParams:
    def test_unit_circle_ellipse(self):
        assert Ellipse(ORIGIN, 1, 1) == UNIT_CIRCLE

    def test_overspecified(self):
        with pytest.raises(ValueError, match="not both"):
            Ellipse(ORIGIN, 1, 1, r1_angle=0, r1_direction=(1, 0))

    def test_wrong_r1_direction(self):
        with pytest.raises(ValueError, match="must be a 2d vector or an ideal point"):
            Ellipse(ORIGIN, 1, 1, r1_direction=(1, 1, 1))

    def test_center(self):
        center = symbols("x,y")
        r1, r2 = symbols("r1,r2")
        r1_direction = symbols("r1_x,r1_y")
        ellipse = Ellipse(center, r1, r2, r1_direction=r1_direction)
        assert center == factor(ConicCenter(ellipse))

    def test_axis_direction_vector_length_invariance(self):
        center = symbols("x,y")
        r1, r2, r1_dir_x, r1_dir_y = symbols("r1,r2,dx,dy", positive=True)
        conic1 = Ellipse(center, r1, r2, r1_direction=(r1_dir_x, r1_dir_y))
        conic2 = Ellipse(center, r1, r2, r1_direction=(r1_dir_x * -2, r1_dir_y * -2))
        assert IsNonZeroMultiple(conic1, conic2)

    def test_axis_lengths_circle(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = Ellipse(center, r, r)
        assert SemiAxisLengths(circle) == (r, r)

    def test_axis_lengths_numeric(self):
        ellipse = Ellipse((1, 2), 3, 4, r1_direction=(5, 6))
        assert SemiAxisLengths(ellipse) == (3, 4)
        ellipse = Ellipse((1, 2), 3, 4, r1_angle=pi / 6)
        assert SemiAxisLengths(ellipse) == (3, 4)

    def test_axis_lengths_general_case(self):
        center = symbols("x,y")
        r_min, r_diff = symbols("r_min,r_diff", positive=True)
        ellipse = Ellipse(center, r_min, r_min + r_diff, r1_direction=(73, -25))
        axes = tuple(factor(simplify(length)) for length in SemiAxisLengths(ellipse))
        assert (r_min, r_min + r_diff) == axes or (r_min + r_diff, r_min) == axes


class TestEllipseFromFociAndPoint:
    def test_numeric_ellipse(self):
        f1 = (1, 2)
        f2 = (3, 4)
        p = (0, 0)
        ellipse = EllipseFromFociAndPoint(f1, f2, p)
        assert IsEllipse(ellipse)
        center = [coord.simplify() for coord in ConicCenter(ellipse)]
        assert center == [2, 3]
        assert ConicContainsPoint(ellipse, p)

    def test_collinear_points(self):
        f1 = (1, 2)
        f2 = (3, 4)
        conic = EllipseFromFociAndPoint(f1, f2, (2, 3))
        assert IsNonZeroMultiple(conic, DoubleLine(LineBetween(f1, f2)))
        conic = EllipseFromFociAndPoint(f1, f2, (1, 2))
        assert IsNonZeroMultiple(conic, DoubleLine(LineBetween(f1, f2)))

    def test_coincident_foci(self):
        f = (2, 3)
        p = (5, 7)
        circle = EllipseFromFociAndPoint(f, f, p)
        assert IsNonZeroMultiple(circle, Circle(f, 5))

    def test_coincident_foci_and_point(self):
        assert EllipseFromFociAndPoint((1, 2), (1, 2), (1, 2)).is_zero_matrix

    def test_ideal_point_focus(self):
        f1 = (1, 2)
        f2 = (3, 4, 0)
        p = (5, 6)
        assert nan in EllipseFromFociAndPoint(f1, f2, p)

    def test_ideal_incident_point(self):
        f1 = (1, 2)
        f2 = (3, 4)
        p = (5, 6, 0)
        assert nan in EllipseFromFociAndPoint(f1, f2, p)


class TestSteinerEllipse:
    def test_circumellipse(self):
        p1 = (4, 1)
        p2 = (7, 3)
        p3 = (5, 5)
        ellipse = SteinerEllipse(p1, p2, p3)
        assert ConicContainsPoint(ellipse, PointToVec3(p1))
        assert ConicContainsPoint(ellipse, PointToVec3(p2))
        assert ConicContainsPoint(ellipse, PointToVec3(p3))
        assert ConicCenter(ellipse) == Centroid(p1, p2, p3)

    def test_inellipse(self):
        p1 = (4, 1)
        p2 = (7, 3)
        p3 = (5, 5)
        ellipse = SteinerInellipse(p1, p2, p3)
        assert ConicContainsPoint(ellipse, PointToVec3(Centroid(p2, p3)))
        assert ConicContainsPoint(ellipse, PointToVec3(Centroid(p3, p1)))
        assert ConicContainsPoint(ellipse, PointToVec3(Centroid(p1, p2)))
        assert ConicCenter(ellipse) == Centroid(p1, p2, p3)
