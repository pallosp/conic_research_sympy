import pytest
from sympy import Matrix, Rational, symbols

from lib.line import IDEAL_LINE
from lib.point import Centroid, IdealPointOnLine


class TestIdealPointOnLine:
    def test_symbolic(self):
        line = Matrix(symbols("a b c"))
        point = IdealPointOnLine(line)
        assert point[2] == 0  # point should be an ideal point
        assert line.dot(point) == 0  # point should be on the line

    def test_ideal_line(self):
        assert IdealPointOnLine(IDEAL_LINE) == Matrix([0, 0, 0])


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
