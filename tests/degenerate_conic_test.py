from sympy import I, Matrix

from lib.circle import Circle
from lib.conic import IdealPoints
from lib.conic_classification import IsDegenerate, IsFiniteConic
from lib.degenerate_conic import LinePair, PointConic, SplitToLines
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, VerticalLine
from lib.matrix import ConicMatrix, IsNonZeroMultiple, IsRealMatrix
from tests.util import AreProjectiveSetsEqual


class TestLinePair:
    def test_line_pair(self):
        double_ideal = ConicMatrix(0, 0, 0, 0, 0, 1)
        assert IsNonZeroMultiple(LinePair(IDEAL_LINE, IDEAL_LINE), double_ideal)
        plus = ConicMatrix(0, 1, 0, 0, 0, 0)
        assert IsNonZeroMultiple(LinePair(X_AXIS, Y_AXIS), plus)


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
