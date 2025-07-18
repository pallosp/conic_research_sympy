from sympy import I, Matrix, pi, symbols
from sympy.abc import x, y

from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import ConicFromPoly
from lib.conic_classification import (
    IsCircular,
    IsComplexEllipse,
    IsDegenerate,
    IsFiniteConic,
    IsHyperbola,
    IsParabola,
)
from lib.degenerate_conic import LinePair
from lib.ellipse import Ellipse
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine
from lib.matrix import ConicMatrix


class TestIsDegenerate:
    def test_numeric(self):
        assert IsDegenerate(Circle((1, 2), 0))
        assert IsDegenerate(UNIT_CIRCLE) is False

    def test_symbolic(self):
        assert IsDegenerate(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None
        assert IsDegenerate(Matrix.diag(symbols("a,c,f", positive=True))) is False
        line_pair = LinePair(Matrix(symbols("a,b,c")), Matrix(symbols("d,e,f")))
        assert IsDegenerate(line_pair)
        point = Circle(symbols("x,y"), 0)
        assert IsDegenerate(point)


class TestIsFiniteConic:
    def test_numeric_point(self):
        assert IsFiniteConic(Circle((1, 2), 0))

    def test_numeric_circle(self):
        assert IsFiniteConic(UNIT_CIRCLE)
        assert IsFiniteConic(Circle((1, 2), 3))

    def test_numeric_parabola(self):
        assert IsFiniteConic(ConicFromPoly(x * x - y)) is False

    def test_numeric_hyperbola(self):
        assert IsFiniteConic(ConicFromPoly(x * y - 1)) is False

    def test_numeric_line_pair(self):
        assert IsFiniteConic(LinePair(X_AXIS, Y_AXIS)) is False

    def test_symbolic_conic(self):
        assert IsFiniteConic(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_circle(self):
        circle = Circle(symbols("x,y"), symbols("r", positive=True))
        assert IsFiniteConic(circle)

    def test_symbolic_complex_ellipse(self):
        complex_ellipse = Matrix.diag(symbols("a,c,f", positive=True))
        assert IsFiniteConic(complex_ellipse)

    def test_symbolic_point(self):
        point = Circle(symbols("x,y"), 0)
        assert IsFiniteConic(point)

    def test_symbolic_line_pair(self):
        line1 = Matrix(symbols("a,b,c", real=True))
        line2 = Matrix(symbols("d,e,f", real=True))
        assert IsFiniteConic(LinePair(line1, line2)) is False


class TestIsComplexEllipse:
    def test_numeric_ellipse(self):
        assert IsComplexEllipse(Ellipse((1, 2), 3, 4)) is False
        assert IsComplexEllipse(Ellipse((1, 2), 3 * I, 4 * I)) is True

    def test_symbolic_ellipse(self):
        center = symbols("x,y")
        r = symbols("r1,r2", positive=True)
        imag_r = [r[0] * I, r[1] * I]
        dir = symbols("dx,dy", positive=True)
        assert IsComplexEllipse(Ellipse(center, *r)) is False
        assert IsComplexEllipse(Ellipse(center, *r, r1_direction=dir)) is False
        assert IsComplexEllipse(Ellipse(center, *imag_r)) is True
        assert IsComplexEllipse(Ellipse(center, *imag_r, r1_direction=dir)) is True

    def test_symbolic_point(self):
        point = Circle(symbols("x,y"), 0)
        assert IsComplexEllipse(point) is False


class TestIsParabola:
    def test_numeric(self):
        assert IsParabola(ConicFromPoly(x * x - y))
        assert not IsParabola(LinePair(X_AXIS, HorizontalLine(1)))
        assert not IsParabola(UNIT_CIRCLE)

    def test_symbolic(self):
        assert IsParabola(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None


class TestIsHyperbola:
    def test_numeric(self):
        assert IsHyperbola(ConicFromPoly(x * y - 1))
        assert IsHyperbola(ConicFromPoly(x * x - y * y - 1))
        assert not IsHyperbola(ConicFromPoly(x * x - y * y))
        assert not IsHyperbola(UNIT_CIRCLE)

    def test_symbolic(self):
        assert IsHyperbola(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None


class TestIsCircular:
    def test_circle(self):
        assert IsCircular(UNIT_CIRCLE)
        assert IsCircular(Circle((1, 2), 3))

    def test_point(self):
        assert IsCircular(Circle((1, 2), 0))

    def test_complex_circle(self):
        assert IsCircular(Circle((1, 2), I))

    def test_ellipse(self):
        assert IsCircular(Ellipse((0, 0), 1, 1))
        assert IsCircular(Ellipse((0, 0), 1, 1, r1_angle=pi / 4))
        assert not IsCircular(Ellipse((0, 0), 1, 2))

    def test_hyperbola(self):
        assert not IsCircular(ConicFromPoly(x * y - 1))

    def test_lines(self):
        assert not IsCircular(LinePair(IDEAL_LINE, IDEAL_LINE))
        assert not IsCircular(LinePair(IDEAL_LINE, X_AXIS))
        assert not IsCircular(LinePair(X_AXIS, X_AXIS))
        assert not IsCircular(LinePair(X_AXIS, Y_AXIS))

    def test_symbolic_conic(self):
        r = symbols("r")
        assert IsCircular(Circle((0, 0), r))
        pos1, pos2 = symbols("pos1,pos2", positive=True)
        neg = symbols("neg", negative=True)
        assert IsCircular(ConicFromPoly(pos1 * x**2 + pos1 * y**2)) is True
        assert IsCircular(ConicFromPoly(pos1 * x**2 + neg * y**2)) is False
        assert IsCircular(ConicFromPoly(pos1 * x**2 + pos2 * y**2)) is None
