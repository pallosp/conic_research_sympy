from sympy import nan

from lib.distance import PointLineDistance
from lib.line import IDEAL_LINE, VerticalLine, X_AXIS
from lib.point import IdealPoint


class TestPointLineDistance:
    def test_euclidean_point(self):
        assert PointLineDistance((1, 2), VerticalLine(4)) == 3
        assert PointLineDistance((5, 6), VerticalLine(4)) == -1

    def test_ideal_point(self):
        assert PointLineDistance(IdealPoint(1, 2), X_AXIS) == nan

    def test_ideal_line(self):
        assert PointLineDistance((1, 2), IDEAL_LINE).is_infinite
