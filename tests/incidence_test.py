from sympy import I, Matrix, symbols
from sympy.abc import x, y

from lib.circle import circle
from lib.conic import conic_from_focus_and_directrix, conic_from_poly
from lib.degenerate_conic import line_pair_conic
from lib.incidence import (
    are_collinear,
    are_concurrent,
    conic_contains_line,
    conic_contains_point,
    line_contains_point,
)
from lib.line import (
    IDEAL_LINE,
    X_AXIS,
    Y_AXIS,
    horizontal_line,
    line_through_point,
)


class TestConicContainsPoint:
    def test_numeric(self):
        conic = conic_from_poly(x * y - 6)
        assert conic_contains_point(conic, (3, 2)) is True
        assert conic_contains_point(conic, (1, 0, 0)) is True
        assert conic_contains_point(conic, (1, 1, 0)) is False

    def test_symbolic(self):
        x, y = symbols("x y", real=True)
        conic = conic_from_poly(x * y - 6, x=x, y=y)
        assert conic_contains_point(conic, (x, y)) is None
        assert conic_contains_point(conic, (x, 0, 0)) is True
        assert conic_contains_point(conic, (x, -x)) is False


class TestConicContainsLine:
    def test_circle_symbolic(self):
        a_circle = circle((1, 2), 3)
        line = Matrix(symbols("a b c", positive=True))
        assert conic_contains_line(a_circle, line) is False

    def test_parabola_symbolic(self):
        focus = (0, 0)
        directrix = Matrix(symbols("d1,d2,d3", positive=True))
        conic = conic_from_focus_and_directrix(focus, directrix, 1)
        line = Matrix(symbols("a,b,c", positive=True))
        assert conic_contains_line(conic, line) is False

    def test_line_pair_symbolic(self):
        line1 = Matrix(symbols("a b c"))
        line2 = X_AXIS
        conic = line_pair_conic(line1, line2)
        assert conic_contains_line(conic, line1) is True
        assert conic_contains_line(conic, line2) is True
        assert conic_contains_line(conic, Y_AXIS) is None
        assert conic_contains_line(conic, IDEAL_LINE) is None

    def test_coincident_line_pair_numeric(self):
        conic = line_pair_conic(X_AXIS, X_AXIS)
        assert conic_contains_line(conic, X_AXIS) is True
        assert conic_contains_line(conic, Y_AXIS) is False
        assert conic_contains_line(conic, horizontal_line(1)) is False

    def test_point_conic_numeric(self):
        conic = circle((0, 0), 0)
        assert conic_contains_line(conic, X_AXIS) is False
        assert conic_contains_line(conic, horizontal_line(1)) is False
        assert conic_contains_line(conic, Matrix([1, I, 0])) is True


class TestLineContainsPoint:
    def test_numeric(self):
        line = Matrix([1, 1, -3])  # x+y=3
        assert line_contains_point(line, (2, 1)) is True
        assert line_contains_point(line, (2, 2)) is False

    def test_symbolic(self):
        line = Matrix([1, 1, -3])  # x+y=3
        x = symbols("x")
        assert line_contains_point(line, (x, 3 - x)) is True
        assert line_contains_point(line, (x, 2 - x)) is False
        assert line_contains_point(line, (x, 1)) is None


class TestAreCollinear:

    def test_less_than_three_points(self):
        assert are_collinear((1, 2)) is True
        assert are_collinear((1, 2), (3, 4)) is True
        assert are_collinear((1, 2), (3, 4, 0)) is True

    def test_three_points_numeric(self):
        assert are_collinear((1, 2), (3, 4), (5, 6)) is True
        assert are_collinear((1, 2), (3, 4), (5, 7)) is False
        assert are_collinear((1, 2), (3, 4), (1, 1, 0)) is True
        assert are_collinear((1, 2), (3, 4), (1, 2, 0)) is False

    def test_three_points_symbolic(self):
        x = symbols("x")
        assert are_collinear((1, 2), (3, 4), (x, x + 1)) is True
        assert are_collinear((1, 2), (3, 4), (x, x)) is False
        assert are_collinear((1, 2), (3, 4), symbols("x y")) is None

    def test_four_points_numeric(self):
        assert are_collinear((1, 2), (1, 2), (3, 4), (5, 6)) is True
        assert are_collinear((1, 2), (1, 2), (3, 4), (5, 7)) is False
        assert are_collinear((1, 2), (3, 4), (5, 6), (7, 8)) is True
        assert are_collinear((1, 2), (3, 4), (5, 6), (7, 9)) is False


class TestAreConcurrent:

    def test_symbolic_lines(self):
        p = symbols("x y")
        x1, y1, x2, y2, x3, y3, x4, y4 = symbols("x1 y1 x2 y2 x3 y3 x4 y4")
        line1 = line_through_point(p, direction=(x1, y1))
        line2 = line_through_point(p, direction=(x2, y2))
        line3 = line_through_point(p, direction=(x3, y3))
        line4 = line_through_point(p, direction=(x4, y4))
        assert are_concurrent(line1, line2, line3) is True
        assert are_concurrent(line1, line2, line3, line4) is True
        assert are_concurrent(line1, line2, X_AXIS) is None

    def test_symbolic_triangle(self):
        line = Matrix(symbols("a b c", positive=True))
        assert are_concurrent(X_AXIS, Y_AXIS, line) is False
        assert are_concurrent(X_AXIS, Y_AXIS, line) is False

    def test_symbolic_parallels(self):
        horiz1 = horizontal_line(symbols("y1", positive=True))
        horiz2 = horizontal_line(symbols("y2", negative=True))
        horiz3 = horizontal_line(symbols("y3"))
        assert are_concurrent(X_AXIS, horiz1, horiz2) is True
        assert are_concurrent(X_AXIS, horiz1, horiz2, horiz3) is True
        assert are_concurrent(horiz1, horiz2, Y_AXIS) is False
        assert are_concurrent(horiz1, horiz2, X_AXIS, Y_AXIS) is False
        # This is actually False, but Sympy 1.14 can't prove it.
        assert are_concurrent(horiz1, horiz2, horiz3, Y_AXIS) is None
