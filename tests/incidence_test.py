from sympy import I, Matrix, symbols
from sympy.abc import x, y

from lib.circle import Circle
from lib.conic import ConicFromFocusAndDirectrix, ConicFromPoly
from lib.degenerate_conic import LinePair
from lib.incidence import (
    AreCollinear,
    ConicContainsLine,
    ConicContainsPoint,
    LineContainsPoint,
)
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine


class TestConicContainsPoint:
    def test_numeric(self):
        conic = ConicFromPoly(x * y - 6)
        assert ConicContainsPoint(conic, (3, 2)) is True
        assert ConicContainsPoint(conic, (1, 0, 0)) is True
        assert ConicContainsPoint(conic, (1, 1, 0)) is False

    def test_symbolic(self):
        x, y = symbols("x y", real=True)
        conic = ConicFromPoly(x * y - 6, x=x, y=y)
        assert ConicContainsPoint(conic, (x, y)) is None
        assert ConicContainsPoint(conic, (x, 0, 0)) is True
        assert ConicContainsPoint(conic, (x, -x)) is False


class TestConicContainsLine:
    def test_circle_symbolic(self):
        circle = Circle((1, 2), 3)
        line = Matrix(symbols("a b c", positive=True))
        assert ConicContainsLine(circle, line) is False

    def test_parabola_symbolic(self):
        focus = (0, 0)
        directrix = Matrix(symbols("d1,d2,d3", positive=True))
        conic = ConicFromFocusAndDirectrix(focus, directrix, 1)
        line = Matrix(symbols("a,b,c", positive=True))
        assert ConicContainsLine(conic, line) is False

    def test_line_pair_symbolic(self):
        line1 = Matrix(symbols("a b c"))
        line2 = X_AXIS
        conic = LinePair(line1, line2)
        assert ConicContainsLine(conic, line1) is True
        assert ConicContainsLine(conic, line2) is True
        assert ConicContainsLine(conic, Y_AXIS) is None
        assert ConicContainsLine(conic, IDEAL_LINE) is None

    def test_coincident_line_pair_numeric(self):
        conic = LinePair(X_AXIS, X_AXIS)
        assert ConicContainsLine(conic, X_AXIS) is True
        assert ConicContainsLine(conic, Y_AXIS) is False
        assert ConicContainsLine(conic, HorizontalLine(1)) is False

    def test_point_conic_numeric(self):
        conic = Circle((0, 0), 0)
        assert ConicContainsLine(conic, X_AXIS) is False
        assert ConicContainsLine(conic, HorizontalLine(1)) is False
        assert ConicContainsLine(conic, Matrix([1, I, 0])) is True


class TestLineContainsPoint:
    def test_numeric(self):
        line = Matrix([1, 1, -3])  # x+y=3
        assert LineContainsPoint(line, (2, 1)) is True
        assert LineContainsPoint(line, (2, 2)) is False

    def test_symbolic(self):
        line = Matrix([1, 1, -3])  # x+y=3
        x = symbols("x")
        assert LineContainsPoint(line, (x, 3 - x)) is True
        assert LineContainsPoint(line, (x, 2 - x)) is False
        assert LineContainsPoint(line, (x, 1)) is None


class TestAreCollinear:

    def test_less_than_three_points(self):
        assert AreCollinear((1, 2)) is True
        assert AreCollinear((1, 2), (3, 4)) is True
        assert AreCollinear((1, 2), (3, 4, 0)) is True

    def test_three_points_numeric(self):
        assert AreCollinear((1, 2), (3, 4), (5, 6)) is True
        assert AreCollinear((1, 2), (3, 4), (5, 7)) is False
        assert AreCollinear((1, 2), (3, 4), (1, 1, 0)) is True
        assert AreCollinear((1, 2), (3, 4), (1, 2, 0)) is False

    def test_three_points_symbolic(self):
        x = symbols("x")
        assert AreCollinear((1, 2), (3, 4), (x, x + 1)) is True
        assert AreCollinear((1, 2), (3, 4), (x, x)) is False
        assert AreCollinear((1, 2), (3, 4), symbols("x y")) is None

    def test_four_points_numeric(self):
        assert AreCollinear((1, 2), (1, 2), (3, 4), (5, 6)) is True
        assert AreCollinear((1, 2), (1, 2), (3, 4), (5, 7)) is False
        assert AreCollinear((1, 2), (3, 4), (5, 6), (7, 8)) is True
        assert AreCollinear((1, 2), (3, 4), (5, 6), (7, 9)) is False
