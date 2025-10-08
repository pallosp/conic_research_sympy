import pytest
from sympy import I, Matrix, symbols

from lib.circle import Circle
from lib.conic import ConicFromPoly, IdealPoints
from lib.conic_classification import IsDegenerate, IsFiniteConic
from lib.degenerate_conic import (
    DoubleLine,
    ExtractPoint,
    LinePair,
    PointConic,
    SplitToLines,
)
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, VerticalLine
from lib.matrix import ConicMatrix, IsNonZeroMultiple, IsRealMatrix
from lib.point import PointToVec3
from tests.utils import AreProjectiveSetsEqual


class TestLinePair:
    def test_invalid_argument(self):
        with pytest.raises(ValueError, match="must be 3-dimensional column vectors"):
            LinePair(Matrix([1, 2, 3]), Matrix([1, 2, 3]).T)
        with pytest.raises(ValueError, match="must be 3-dimensional column vectors"):
            LinePair(Matrix([1, 2, 3]).T, Matrix([1, 2, 3]))

    def test_numeric_line_pair(self):
        double_ideal = ConicMatrix(0, 0, 0, 0, 0, 1)
        assert IsNonZeroMultiple(LinePair(IDEAL_LINE, IDEAL_LINE), double_ideal)
        plus = ConicMatrix(0, 1, 0, 0, 0, 0)
        assert IsNonZeroMultiple(LinePair(X_AXIS, Y_AXIS), plus)


class TestDoubleLine:
    def test_invalid_argument(self):
        with pytest.raises(ValueError, match="must be a 3-dimensional column vector"):
            DoubleLine(Matrix([1, 2, 3]).T)

    def test_numeric_line_pair(self):
        assert DoubleLine(X_AXIS) == ConicFromPoly(symbols("y") ** 2)


class TestPointConic:
    def test_real_point(self):
        point = (1, 2)
        conic = PointConic(point)
        assert IsDegenerate(conic) is True
        assert IsFiniteConic(conic) is True

    def test_ideal_point(self):
        point = (1, 2, 0)
        conic = PointConic(point)
        assert IsDegenerate(conic) is True
        assert AreProjectiveSetsEqual(IdealPoints(conic), [point, point])
        assert IsRealMatrix(SplitToLines(conic)[0]) is False


class TestSplitToLines:
    def test_horizontal_lines(self):
        lines = [HorizontalLine(1), HorizontalLine(2)]
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_vertical_lines(self):
        lines = [VerticalLine(1), VerticalLine(2)]
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_double_line(self):
        lines = [Matrix([1, 2, 3])] * 2
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_double_ideal_line(self):
        lines = [IDEAL_LINE] * 2
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_general_case(self):
        lines = [Matrix([1, 2, 3]), Matrix([4, 5, 6])]
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_point_conic(self):
        conic = Circle((0, 0), 0)
        assert AreProjectiveSetsEqual(
            [Matrix([1, I, 0]), Matrix([1, -I, 0])],
            SplitToLines(conic),
        )


class TestExtractPoint:
    def test_symbolic_real_point_conic(self):
        point = symbols("x,y", real=True)
        zero_circle = Circle(point, 0)
        assert IsNonZeroMultiple(ExtractPoint(zero_circle), PointToVec3(point))

        finite_point_conic = PointConic(point)
        assert IsNonZeroMultiple(ExtractPoint(finite_point_conic), PointToVec3(point))

    def test_symbolic_ideal_point_conic(self):
        point = [*symbols("x,y", positive=True), 0]
        conic = PointConic(point)
        assert IsNonZeroMultiple(ExtractPoint(conic), PointToVec3(point))

    def test_zero_matrix(self):
        conic = Matrix.zeros(3, 3)
        assert ExtractPoint(conic) == Matrix([0, 0, 0])

    def test_numeric_line_pairs(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        conic = LinePair(line1, line2)
        intersection = line1.cross(line2)
        assert IsNonZeroMultiple(ExtractPoint(conic), intersection)

    def test_symbolic_double_line(self):
        line = Matrix(symbols("a b c"))
        conic = LinePair(line, line)
        assert ExtractPoint(conic) == Matrix([0, 0, 0])
