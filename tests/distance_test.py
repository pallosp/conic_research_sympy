import pytest
from sympy import Matrix, Rational, nan, symbols

from lib.distance import (
    parallel_line_distance,
    point_line_distance,
    point_point_distance,
)
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, horizontal_line, vertical_line
from lib.point import ideal_point


class TestPointPointDistance:
    def test_two_euclidean_points(self):
        assert point_point_distance((2, 3), (5, 7)) == 5
        assert point_point_distance((2, 3, 1), (5, 7, 1)) == 5
        assert point_point_distance((-2, -3, -1), (5, 7, 1)) == -5

    def test_one_ideal_point(self):
        assert point_point_distance((1, 2), ideal_point(3, 4)).is_infinite

    def test_two_ideal_points(self):
        assert point_point_distance(ideal_point(1, 2), ideal_point(3, 4)) == nan
        assert point_point_distance(ideal_point(1, 2), ideal_point(1, 2)) == nan


class TestPointLineDistance:
    def test_euclidean_point(self):
        assert point_line_distance((1, 2), vertical_line(4)) == 3
        assert point_line_distance((5, 6), vertical_line(4)) == -1

    def test_ideal_point_real_line(self):
        assert point_line_distance(ideal_point(1, 2), X_AXIS) == nan

    def test_real_point_ideal_line(self):
        assert point_line_distance((1, 2), IDEAL_LINE).is_infinite

    def test_ideal_point_ideal_line(self):
        assert point_line_distance(ideal_point(1, 2), IDEAL_LINE) == nan


class TestParallelLineDistance:
    def test_parallel_horizontal_lines(self):
        line1 = horizontal_line(2)
        line2 = horizontal_line(5)
        assert parallel_line_distance(line1, line2) == 3

    def test_parallel_vertical_lines(self):
        line1 = vertical_line(8)
        line2 = vertical_line(21)
        assert parallel_line_distance(line1, line2) == 13

    def test_parallel_lines(self):
        line1 = Matrix([3, 4, 5])
        line2 = Matrix([3, 4, 6])
        assert parallel_line_distance(line1, line2) == Rational(1, 5)
        assert parallel_line_distance(line2, line1) == Rational(1, 5)
        assert parallel_line_distance(line1, -line2) == -Rational(1, 5)

    def test_crossing_lines(self):
        with pytest.raises(ValueError, match="The lines must be parallel"):
            parallel_line_distance(X_AXIS, Y_AXIS)

    def test_one_ideal_line(self):
        assert parallel_line_distance(X_AXIS, IDEAL_LINE).is_infinite

    def test_two_ideal_lines(self):
        assert parallel_line_distance(IDEAL_LINE, IDEAL_LINE) == nan

    def test_symbolic_one_ideal_line(self):
        finite_line = Matrix(symbols("a b c", positive=True))
        assert parallel_line_distance(finite_line, IDEAL_LINE).is_infinite

    def test_symbolic_crossing_lines(self):
        non_horizontal_line = Matrix(symbols("a b c", positive=True))
        with pytest.raises(ValueError, match="The lines must be parallel"):
            parallel_line_distance(X_AXIS, non_horizontal_line)

    def test_symbolic_undecidable_parallelness(self):
        line = Matrix(symbols("a b c"))
        parallel_line_distance(line, X_AXIS)  # shouldn't raise a ValueError
