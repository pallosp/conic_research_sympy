from lib.degenerate_conic import LinePair
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS
from lib.matrix import ConicMatrix, IsNonZeroMultiple


class TestLinePair:
    def test_line_pair(self):
        double_ideal = ConicMatrix(0, 0, 0, 0, 0, 1)
        assert IsNonZeroMultiple(LinePair(IDEAL_LINE, IDEAL_LINE), double_ideal)
        plus = ConicMatrix(0, 1, 0, 0, 0, 0)
        assert IsNonZeroMultiple(LinePair(X_AXIS, Y_AXIS), plus)
