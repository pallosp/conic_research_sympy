from sympy import Function, I, Matrix, nan, sqrt, symbols

from lib.circle import UNIT_CIRCLE
from lib.degenerate_conic import LinePair
from lib.intersection import ConicXLine, LineXLine
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, VerticalLine
from lib.matrix import ConicMatrix, IsNonZeroMultiple, QuadraticForm
from lib.point import ORIGIN, IdealPoint, PointToXY


class TestLineXLine:
    def test_horizonal_x_vertical(self):
        x, y = symbols("x,y")
        line1 = VerticalLine(x)
        line2 = HorizontalLine(y)
        assert IsNonZeroMultiple(LineXLine(line1, line2), [x, y, 1])

    def test_parallel_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        assert LineXLine(line1, line2) == Matrix([2, -1, 0])

    def test_coincident_lines(self):
        line = Matrix([1, 2, 3])
        assert LineXLine(line, line) == Matrix([0, 0, 0])


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
        assert IsNonZeroMultiple(p1, IdealPoint(0, 1))
        assert IsNonZeroMultiple(p2, IdealPoint(1, 0))

    def test_x_axis_plus_ideal_line_conic(self):
        conic = LinePair(X_AXIS, IDEAL_LINE)
        p1, p2 = ConicXLine(conic, Y_AXIS)
        assert IsNonZeroMultiple(p1, IdealPoint(0, 1))
        assert IsNonZeroMultiple(p2, ORIGIN)

    def test_y_axis_plus_ideal_line_conic(self):
        conic = LinePair(Y_AXIS, IDEAL_LINE)
        p1, p2 = ConicXLine(conic, X_AXIS)
        assert IsNonZeroMultiple(p1, IdealPoint(1, 0))
        assert IsNonZeroMultiple(p2, ORIGIN)

    def test_double_x_axis_conic(self):
        conic = LinePair(X_AXIS, X_AXIS)
        assert IsNonZeroMultiple(ConicXLine(conic, Y_AXIS)[0], Matrix([0, 0, 1]))
        assert IsNonZeroMultiple(ConicXLine(conic, IDEAL_LINE)[0], Matrix([1, 0, 0]))

    def test_conic_contains_line(self):
        assert ConicXLine(LinePair(X_AXIS, Y_AXIS), X_AXIS) is nan
        assert ConicXLine(LinePair(X_AXIS, X_AXIS), X_AXIS) is nan

    def test_zero_line(self):
        zero_line = Matrix([0, 0, 0])
        assert ConicXLine(LinePair(X_AXIS, Y_AXIS), zero_line) is nan

    def test_zero_conic(self):
        zero_conic = Matrix.zeros(3, 3)
        assert ConicXLine(zero_conic, Matrix([1, 2, 3])) is nan

    def test_double_intersection(self):
        intersections = ConicXLine(UNIT_CIRCLE, HorizontalLine(1))
        assert IsNonZeroMultiple(intersections[0], Matrix([0, 1, 1]))
        assert IsNonZeroMultiple(intersections[1], Matrix([0, 1, 1]))

    def test_complex_intersection(self):
        intersections = ConicXLine(UNIT_CIRCLE, HorizontalLine(2))
        intersections = sorted((PointToXY(i) for i in intersections), key=str)
        assert intersections == sorted([(-sqrt(3) * I, 2), (sqrt(3) * I, 2)], key=str)

    def test_symbolic_conic(self):
        conic = ConicMatrix(*symbols("a b c d e f"))
        intersections = ConicXLine(conic, IDEAL_LINE)
        assert isinstance(intersections, Function)
