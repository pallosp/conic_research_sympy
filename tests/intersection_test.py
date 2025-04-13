from sympy import Matrix, nan

from lib.degenerate_conic import LinePair
from lib.intersection import ConicXLine
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS
from lib.matrix import ConicMatrix, IsScalarMultiple, QuadraticForm


class TestConicXLine:
    def test_general_case(self):
        conic = ConicMatrix(1, 2, 3, 4, 5, 6)
        line = Matrix([1, 2, 3])
        p1, p2 = ConicXLine(conic, line)
        assert QuadraticForm(conic, p1).equals(0)
        assert QuadraticForm(conic, p2).equals(0)

    def test_x_y_axes_conic(self):
        conic = LinePair(X_AXIS, Y_AXIS)
        p1, p2 = ConicXLine(conic, IDEAL_LINE)
        assert IsScalarMultiple(p1, Matrix([1, 0, 0]))
        assert IsScalarMultiple(p2, Matrix([0, 1, 0]))

    def test_x_axis_plus_ideal_line_conic(self):
        conic = LinePair(X_AXIS, IDEAL_LINE)
        p1, p2 = ConicXLine(conic, Y_AXIS)
        assert IsScalarMultiple(p1, Matrix([0, 1, 0]))
        assert IsScalarMultiple(p2, Matrix([0, 0, 1]))

    def test_y_axis_plus_ideal_line_conic(self):
        conic = LinePair(Y_AXIS, IDEAL_LINE)
        p1, p2 = ConicXLine(conic, X_AXIS)
        assert IsScalarMultiple(p1, Matrix([0, 0, 1]))
        assert IsScalarMultiple(p2, Matrix([1, 0, 0]))

    def test_double_x_axis_conic(self):
        conic = LinePair(X_AXIS, X_AXIS)
        assert IsScalarMultiple(ConicXLine(conic, Y_AXIS)[0], Matrix([0, 0, 1]))
        assert IsScalarMultiple(ConicXLine(conic, IDEAL_LINE)[0], Matrix([1, 0, 0]))

    def test_conic_containing_line(self):
        assert ConicXLine(LinePair(X_AXIS, Y_AXIS), X_AXIS) == nan
        assert ConicXLine(LinePair(X_AXIS, X_AXIS), X_AXIS) == nan

    def test_zero_line(self):
        zero_line = Matrix([0, 0, 0])
        assert ConicXLine(LinePair(X_AXIS, Y_AXIS), zero_line) == nan

    def test_zero_conic(self):
        zero_conic = Matrix.zeros(3, 3)
        assert ConicXLine(zero_conic, Matrix([1, 2, 3])) == nan
