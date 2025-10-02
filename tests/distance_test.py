import pytest
from sympy import Matrix, Rational, nan, symbols

from lib.distance import ParallelLineDistance, PointLineDistance, PointPointDistance
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, VerticalLine
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


class TestParallelLineDistance:
    def test_parallel_horizontal_lines(self):
        line1 = HorizontalLine(2)
        line2 = HorizontalLine(5)
        assert ParallelLineDistance(line1, line2) == 3

    def test_parallel_vertical_lines(self):
        line1 = VerticalLine(8)
        line2 = VerticalLine(21)
        assert ParallelLineDistance(line1, line2) == 13

    def test_parallel_lines(self):
        line1 = Matrix([3, 4, 5])
        line2 = Matrix([3, 4, 6])
        assert ParallelLineDistance(line1, line2) == Rational(1, 5)
        assert ParallelLineDistance(line2, line1) == Rational(1, 5)
        assert ParallelLineDistance(line1, -line2) == -Rational(1, 5)

    def test_crossing_lines(self):
        with pytest.raises(ValueError, match="The lines must be parallel"):
            ParallelLineDistance(X_AXIS, Y_AXIS)

    def test_one_ideal_line(self):
        assert ParallelLineDistance(X_AXIS, IDEAL_LINE).is_infinite

    def test_two_ideal_lines(self):
        assert ParallelLineDistance(IDEAL_LINE, IDEAL_LINE) == nan

    def test_symbolic_one_ideal_line(self):
        finite_line = Matrix(symbols("a b c", positive=True))
        assert ParallelLineDistance(finite_line, IDEAL_LINE).is_infinite

    def test_symbolic_crossing_lines(self):
        non_horizontal_line = Matrix(symbols("a b c", positive=True))
        with pytest.raises(ValueError, match="The lines must be parallel"):
            ParallelLineDistance(X_AXIS, non_horizontal_line)

    def test_symbolic_undecidable_parallelness(self):
        line = Matrix(symbols("a b c"))
        ParallelLineDistance(line, X_AXIS)  # shouldn't raise a ValueError
