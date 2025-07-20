from sympy import I, Matrix, Rational, pi, symbols
from sympy.abc import x, y

from lib.circle import COMPLEX_UNIT_CIRCLE, UNIT_CIRCLE, Circle
from lib.conic import ConicFromFocusAndDirectrix, ConicFromPoly
from lib.conic_classification import (
    IsCircular,
    IsComplexEllipse,
    IsDegenerate,
    IsDoubleLine,
    IsEllipse,
    IsFiniteConic,
    IsHyperbola,
    IsLinePair,
    IsParabola,
)
from lib.degenerate_conic import LinePair
from lib.ellipse import Ellipse
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine
from lib.matrix import ConicMatrix


class TestIsDegenerate:
    def test_numeric(self):
        assert IsDegenerate(Circle((1, 2), 0)) is True
        assert IsDegenerate(UNIT_CIRCLE) is False

    def test_symbolic(self):
        assert IsDegenerate(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None
        assert IsDegenerate(Matrix.diag(symbols("a,c,f", positive=True))) is False
        line_pair = LinePair(Matrix(symbols("a,b,c")), Matrix(symbols("d,e,f")))
        assert IsDegenerate(line_pair) is True
        point = Circle(symbols("x,y"), 0)
        assert IsDegenerate(point) is True


class TestIsFiniteConic:
    def test_numeric_point(self):
        assert IsFiniteConic(Circle((1, 2), 0)) is True

    def test_numeric_circle(self):
        assert IsFiniteConic(UNIT_CIRCLE) is True
        assert IsFiniteConic(Circle((1, 2), 3)) is True

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
        assert IsFiniteConic(circle) is True

    def test_symbolic_complex_ellipse(self):
        complex_ellipse = Matrix.diag(symbols("a,c,f", positive=True))
        assert IsFiniteConic(complex_ellipse) is True

    def test_symbolic_point(self):
        point = Circle(symbols("x,y"), 0)
        assert IsFiniteConic(point) is True

    def test_symbolic_line_pair(self):
        line1 = Matrix(symbols("a,b,c", real=True))
        line2 = Matrix(symbols("d,e,f", real=True))
        assert IsFiniteConic(LinePair(line1, line2)) is False


class TestIsEllipse:
    def test_numeric_circle(self):
        assert IsEllipse(UNIT_CIRCLE) is True
        assert IsEllipse(-UNIT_CIRCLE) is True
        assert IsEllipse(Circle((1, 2), 0)) is False

    def test_numeric_ellipse(self):
        assert IsEllipse(Ellipse((1, 2), 3, 4)) is True

    def test_numeric_hyperbola(self):
        assert IsEllipse(ConicFromPoly(x * y - 1)) is False

    def test_numeric_point(self):
        assert IsEllipse(Circle((1, 2), 0)) is False

    def test_complex_circle(self):
        assert IsEllipse(COMPLEX_UNIT_CIRCLE) is False

    def test_undecidable(self):
        assert IsEllipse(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_ellipse(self):
        center = symbols("x,y")
        assert IsEllipse(Ellipse(center, *symbols("r1,r2", positive=True))) is True
        assert IsEllipse(Ellipse(center, *symbols("r1,r2"))) is None
        assert IsEllipse(COMPLEX_UNIT_CIRCLE) is False
        assert IsEllipse(-COMPLEX_UNIT_CIRCLE) is False


class TestIsComplexEllipse:
    def test_complex_unit_circle(self):
        assert IsComplexEllipse(COMPLEX_UNIT_CIRCLE)
        assert IsComplexEllipse(-COMPLEX_UNIT_CIRCLE)

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
        assert IsParabola(LinePair(X_AXIS, HorizontalLine(1))) is False
        assert IsParabola(UNIT_CIRCLE) is False

    def test_symbolic(self):
        assert IsParabola(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None


class TestIsHyperbola:
    def test_numeric(self):
        assert IsHyperbola(ConicFromPoly(x * y - 1)) is True
        assert IsHyperbola(ConicFromPoly(x * x - y * y - 1)) is True
        assert IsHyperbola(ConicFromPoly(x * x - y * y)) is False
        assert IsHyperbola(UNIT_CIRCLE) is False

    def test_unspecified(self):
        assert IsHyperbola(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_ellipse_of_parabola(self):
        focus = (0, 0)
        # Real line not going through the focus
        directrix = Matrix(symbols("a,b,c", positive=True))

        circle = ConicFromFocusAndDirectrix(focus, directrix, 0)
        assert IsHyperbola(circle) is False
        ellipse = ConicFromFocusAndDirectrix(focus, directrix, Rational(1, 2))
        assert IsHyperbola(ellipse) is False
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert IsHyperbola(parabola) is False

    def test_symbolic_hyperbola(self):
        focus = (0, 0)
        # Real line not going through the focus
        directrix = Matrix(symbols("a,b,c", positive=True))
        # Eccentricity greater than 1
        eccentricity = symbols("e", positive=True) + 1
        hyperbola = ConicFromFocusAndDirectrix(focus, directrix, eccentricity)
        assert IsHyperbola(hyperbola) is True


class TestIsCircular:
    def test_circle(self):
        assert IsCircular(UNIT_CIRCLE) is True
        assert IsCircular(Circle((1, 2), 3)) is True

    def test_point(self):
        assert IsCircular(Circle((1, 2), 0)) is True

    def test_complex_circle(self):
        assert IsCircular(Circle((1, 2), I)) is True

    def test_ellipse(self):
        assert IsCircular(Ellipse((0, 0), 1, 1)) is True
        assert IsCircular(Ellipse((0, 0), 1, 1, r1_angle=pi / 4)) is True
        assert IsCircular(Ellipse((0, 0), 1, 2)) is False

    def test_hyperbola(self):
        assert IsCircular(ConicFromPoly(x * y - 1)) is False

    def test_lines(self):
        assert IsCircular(LinePair(IDEAL_LINE, IDEAL_LINE)) is False
        assert IsCircular(LinePair(IDEAL_LINE, X_AXIS)) is False
        assert IsCircular(LinePair(X_AXIS, X_AXIS)) is False
        assert IsCircular(LinePair(X_AXIS, Y_AXIS)) is False

    def test_symbolic_conic(self):
        r = symbols("r")
        assert IsCircular(Circle((0, 0), r)) is True
        pos1, pos2 = symbols("pos1,pos2", positive=True)
        neg = symbols("neg", negative=True)
        assert IsCircular(ConicFromPoly(pos1 * x**2 + pos1 * y**2)) is True
        assert IsCircular(ConicFromPoly(pos1 * x**2 + neg * y**2)) is False
        assert IsCircular(ConicFromPoly(pos1 * x**2 + pos2 * y**2)) is None


class TestIsLinePair:
    def test_zero_matrix(self):
        assert IsLinePair(Matrix.zeros(3, 3)) is False

    def test_numeric_line_pair(self):
        assert IsLinePair(LinePair(X_AXIS, X_AXIS)) is True
        assert IsLinePair(LinePair(X_AXIS, Y_AXIS)) is True
        assert IsLinePair(LinePair(X_AXIS, HorizontalLine(1))) is True
        assert IsLinePair(LinePair(X_AXIS, IDEAL_LINE)) is True
        assert IsLinePair(LinePair(IDEAL_LINE, IDEAL_LINE)) is True
        assert IsLinePair(ConicFromPoly(x * x - y * y)) is True

    def test_numeric_hyperbola(self):
        assert IsLinePair(ConicFromPoly(x * y - 1)) is False

    def test_numeric_parabola(self):
        assert IsLinePair(ConicFromPoly(x * x - y)) is False

    def test_numeric_circle(self):
        assert IsLinePair(UNIT_CIRCLE) is False

    def test_undecidable(self):
        assert IsLinePair(ConicMatrix(*symbols("a,b,c,d,e,f", positive=True))) is None

    def test_symbolic_line_pair(self):
        line_or_zero_vector = Matrix(symbols("a,b,c", real=True))
        line1 = Matrix(symbols("d,e,f", positive=True))
        line2 = Matrix(symbols("g,h,i", positive=True))
        assert IsLinePair(LinePair(line_or_zero_vector, line1)) is None
        assert IsLinePair(LinePair(line1, line1)) is True
        assert IsLinePair(LinePair(line1, line2)) is True

    def test_symbolic_point(self):
        assert IsLinePair(Circle(symbols("x,y"), 0)) is False


class TestIsDoubleLine:
    def test_zero_matrix(self):
        assert IsDoubleLine(Matrix.zeros(3, 3)) is False

    def test_numeric_line_pair(self):
        assert IsDoubleLine(LinePair(X_AXIS, X_AXIS)) is True
        assert IsDoubleLine(LinePair(X_AXIS, Y_AXIS)) is False
        assert IsDoubleLine(LinePair(X_AXIS, HorizontalLine(1))) is False
        assert IsDoubleLine(LinePair(X_AXIS, IDEAL_LINE)) is False
        assert IsDoubleLine(LinePair(IDEAL_LINE, IDEAL_LINE)) is True

    def test_symbolic_line_pair(self):
        line_or_zeros = Matrix(symbols("a,b,c", real=True))
        real_line1 = Matrix(symbols("d,e,f", positive=True))
        real_line2 = Matrix(symbols("g,h,i", positive=True))
        assert IsDoubleLine(LinePair(real_line1, real_line1)) is True
        assert IsDoubleLine(LinePair(line_or_zeros, line_or_zeros)) is None
        assert IsDoubleLine(LinePair(real_line1, real_line2)) is None
        assert IsDoubleLine(LinePair(real_line1, IDEAL_LINE)) is False

    def test_symbolic_point(self):
        assert IsDoubleLine(Circle(symbols("x,y", real=True), 0)) is False
