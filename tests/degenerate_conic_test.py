from sympy import I, Matrix

from lib.circle import Circle
from lib.degenerate_conic import LinePair, SplitToLines
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, VerticalLine
from lib.matrix import ConicMatrix, IsNonZeroMultiple
from tests.util import AreProjectiveSetsEqual


class TestLinePair:
    def test_line_pair(self):
        double_ideal = ConicMatrix(0, 0, 0, 0, 0, 1)
        assert IsNonZeroMultiple(LinePair(IDEAL_LINE, IDEAL_LINE), double_ideal)
        plus = ConicMatrix(0, 1, 0, 0, 0, 0)
        assert IsNonZeroMultiple(LinePair(X_AXIS, Y_AXIS), plus)


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
            [Matrix([1, I, 0]), Matrix([1, -I, 0])], SplitToLines(conic)
        )
