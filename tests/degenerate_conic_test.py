import pytest
from sympy import I, Matrix, symbols

from lib.circle import circle
from lib.conic import IdealPoints, conic_from_poly
from lib.conic_classification import is_degenerate, is_finite_conic
from lib.degenerate_conic import (
    ExtractPoint,
    SplitToLines,
    double_line_conic,
    line_pair_conic,
    point_conic,
)
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, horizontal_line, vertical_line
from lib.matrix import conic_matrix, is_nonzero_multiple, is_real_matrix
from lib.point import point_to_vec3
from tests.utils import are_projective_sets_equal


class TestLinePair:
    def test_invalid_argument(self):
        with pytest.raises(ValueError, match="must be 3-dimensional column vectors"):
            line_pair_conic(Matrix([1, 2, 3]), Matrix([1, 2, 3]).T)
        with pytest.raises(ValueError, match="must be 3-dimensional column vectors"):
            line_pair_conic(Matrix([1, 2, 3]).T, Matrix([1, 2, 3]))

    def test_numeric_line_pair(self):
        double_ideal = conic_matrix(0, 0, 0, 0, 0, 1)
        assert is_nonzero_multiple(
            line_pair_conic(IDEAL_LINE, IDEAL_LINE),
            double_ideal,
        )
        plus = conic_matrix(0, 1, 0, 0, 0, 0)
        assert is_nonzero_multiple(line_pair_conic(X_AXIS, Y_AXIS), plus)


class TestDoubleLine:
    def test_invalid_argument(self):
        with pytest.raises(ValueError, match="must be a 3-dimensional column vector"):
            double_line_conic(Matrix([1, 2, 3]).T)

    def test_numeric_line_pair(self):
        assert double_line_conic(X_AXIS) == conic_from_poly(symbols("y") ** 2)


class TestPointConic:
    def test_real_point(self):
        point = (1, 2)
        conic = point_conic(point)
        assert is_degenerate(conic) is True
        assert is_finite_conic(conic) is True

    def test_ideal_point(self):
        point = (1, 2, 0)
        conic = point_conic(point)
        assert is_degenerate(conic) is True
        assert are_projective_sets_equal(IdealPoints(conic), [point, point])
        assert is_real_matrix(SplitToLines(conic)[0]) is False


class TestSplitToLines:
    def test_horizontal_lines(self):
        lines = [horizontal_line(1), horizontal_line(2)]
        conic = line_pair_conic(*lines)
        assert are_projective_sets_equal(lines, SplitToLines(conic))

    def test_vertical_lines(self):
        lines = [vertical_line(1), vertical_line(2)]
        conic = line_pair_conic(*lines)
        assert are_projective_sets_equal(lines, SplitToLines(conic))

    def test_double_line(self):
        lines = [Matrix([1, 2, 3])] * 2
        conic = line_pair_conic(*lines)
        assert are_projective_sets_equal(lines, SplitToLines(conic))

    def test_double_ideal_line(self):
        lines = [IDEAL_LINE] * 2
        conic = line_pair_conic(*lines)
        assert are_projective_sets_equal(lines, SplitToLines(conic))

    def test_general_case(self):
        lines = [Matrix([1, 2, 3]), Matrix([4, 5, 6])]
        conic = line_pair_conic(*lines)
        assert are_projective_sets_equal(lines, SplitToLines(conic))

    def test_point_conic(self):
        conic = circle((0, 0), 0)
        assert are_projective_sets_equal(
            [Matrix([1, I, 0]), Matrix([1, -I, 0])],
            SplitToLines(conic),
        )


class TestExtractPoint:
    def test_symbolic_real_point_conic(self):
        point = symbols("x,y", real=True)
        zero_circle = circle(point, 0)
        assert is_nonzero_multiple(ExtractPoint(zero_circle), point_to_vec3(point))

        finite_point_conic = point_conic(point)
        assert is_nonzero_multiple(
            ExtractPoint(finite_point_conic),
            point_to_vec3(point),
        )

    def test_symbolic_ideal_point_conic(self):
        point = [*symbols("x,y", positive=True), 0]
        conic = point_conic(point)
        assert is_nonzero_multiple(ExtractPoint(conic), point_to_vec3(point))

    def test_zero_matrix(self):
        conic = Matrix.zeros(3, 3)
        assert ExtractPoint(conic) == Matrix([0, 0, 0])

    def test_numeric_line_pairs(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        conic = line_pair_conic(line1, line2)
        intersection = line1.cross(line2)
        assert is_nonzero_multiple(ExtractPoint(conic), intersection)

    def test_symbolic_double_line(self):
        line = Matrix(symbols("a b c"))
        conic = line_pair_conic(line, line)
        assert ExtractPoint(conic) == Matrix([0, 0, 0])
