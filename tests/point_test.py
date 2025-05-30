import pytest
from sympy import Rational

from lib.point import Centroid


class TestCentroid:
    def test_no_points(self):
        with pytest.raises(AssertionError):
            Centroid()

    def test_one_point(self):
        assert Centroid((1, 2)) == (1, 2)
        assert Centroid((6, 4, 2)) == (3, 2)

    def test_two_points(self):
        assert Centroid((1, 2), (3, 4)) == (2, 3)

    def test_full_precision(self):
        assert Centroid((0, 0), (1, 0), (1, 1)) == (Rational(2, 3), Rational(1, 3))
