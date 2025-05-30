import pytest

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
