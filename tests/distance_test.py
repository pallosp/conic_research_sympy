from sympy import nan

from lib.distance import PointLineDistance, PointPointDistance
from lib.line import IDEAL_LINE, VerticalLine, X_AXIS
from lib.point import IdealPoint


class TestPointPointDistance:
    def test_two_euclidean_points(self):
        assert PointPointDistance((2, 3), (5, 7)) == 5
        assert PointPointDistance((2, 3, 1), (5, 7, 1)) == 5
        assert PointPointDistance((-2, -3, -1), (5, 7, 1)) == -5

    def test_one_ideal_point(self):
        assert PointPointDistance((1, 2), IdealPoint(3, 4)).is_infinite

    def test_two_ideal_points(self):
        assert PointPointDistance(IdealPoint(1, 2), IdealPoint(3, 4)) == nan
        assert PointPointDistance(IdealPoint(1, 2), IdealPoint(1, 2)) == nan


class TestPointLineDistance:
    def test_euclidean_point(self):
        assert PointLineDistance((1, 2), VerticalLine(4)) == 3
        assert PointLineDistance((5, 6), VerticalLine(4)) == -1

    def test_ideal_point_real_line(self):
        assert PointLineDistance(IdealPoint(1, 2), X_AXIS) == nan

    def test_real_point_ideal_line(self):
        assert PointLineDistance((1, 2), IDEAL_LINE).is_infinite

    def test_ideal_point_ideal_line(self):
        assert PointLineDistance(IdealPoint(1, 2), IDEAL_LINE) == nan
