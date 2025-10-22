import pytest
from sympy import Matrix, Rational, nan, symbols

from lib.line import IDEAL_LINE
from lib.matrix import is_nonzero_multiple
from lib.point import centroid, ideal_point_on_line, perpendicular_foot


class TestIdealPointOnLine:
    def test_symbolic(self):
        line = Matrix(symbols("a b c"))
        point = ideal_point_on_line(line)
        assert point[2] == 0  # point should be an ideal point
        assert line.dot(point) == 0  # point should be on the line

    def test_ideal_line(self):
        assert ideal_point_on_line(IDEAL_LINE) == Matrix([0, 0, 0])


class TestCentroid:
    def test_no_points(self):
        with pytest.raises(ValueError, match="At least one point"):
            centroid()

    def test_one_point(self):
        assert centroid((1, 2)) == Matrix([1, 2])
        assert centroid((6, 4, 2)) == Matrix([3, 2])

    def test_two_points(self):
        assert centroid((1, 2), (3, 4)) == Matrix([2, 3])

    def test_full_precision(self):
        expected = Matrix([Rational(2, 3), Rational(1, 3)])
        assert centroid((0, 0), (1, 0), (1, 1)) == expected


class TestPerpendicularFoot:
    def test_finite_point_and_line(self):
        point = (1, -1)
        line = Matrix([1, 1, -4])  # x+y=4
        assert is_nonzero_multiple(perpendicular_foot(point, line), (3, 1))

    def test_ideal_point_or_line(self):
        nan_point = Matrix([nan, nan])
        assert perpendicular_foot((1, 2), IDEAL_LINE) == nan_point
        assert perpendicular_foot((1, 2, 0), Matrix([1, 2, 3])) == nan_point
        assert perpendicular_foot((1, 2, 0), IDEAL_LINE) == nan_point
