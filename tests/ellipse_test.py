from sympy import factor, pi, simplify, symbols

from lib.central_conic import ConicCenter, SemiAxisLengths
from lib.circle import UNIT_CIRCLE
from lib.conic import ConicContainsPoint
from lib.ellipse import Ellipse, SteinerEllipse, SteinerInellipse
from lib.matrix import IsNonZeroMultiple
from lib.point import ORIGIN, Centroid, PointToVec3


class TestEllipseFromParams:
    def test_unit_circle_ellipse(self):
        assert Ellipse(ORIGIN, 1, 1) == UNIT_CIRCLE

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
