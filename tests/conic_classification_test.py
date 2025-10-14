from sympy import Abs, I, Matrix, Rational, pi, symbols
from sympy.abc import x, y

from lib.circle import COMPLEX_UNIT_CIRCLE, UNIT_CIRCLE, Circle
from lib.conic import ConicFromFocusAndDirectrix, ConicFromPoly, FocalAxisDirection
from lib.conic_classification import (
    ConicNormFactor,
    IsCircle,
    IsCircular,
    IsComplexEllipse,
    IsDegenerate,
    IsDoubleLine,
    IsEllipse,
    IsFiniteConic,
    IsFinitePointConic,
    IsHyperbola,
    IsLinePair,
    IsParabola,
    IsPointConic,
)
from lib.degenerate_conic import LinePair, PointConic
from lib.ellipse import Ellipse
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine
from lib.matrix import ConicMatrix, QuadraticForm
from lib.point import ORIGIN, PointToVec3
from lib.transform import TransformConic, Translate


class TestConicNormFactor:
    def test_parabola_with_focus_at_origin(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert ConicNormFactor(parabola) == 1
        assert ConicNormFactor(-parabola) == -1

    def test_general_parabola(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        translation = Translate(*symbols("dx dy", real=True))
        parabola = TransformConic(parabola, translation)
        assert ConicNormFactor(parabola) == 1

    def test_circle(self):
        center = symbols("x y")
        circle = Circle(center, symbols("r", positive=True))
        assert QuadraticForm(circle, PointToVec3(center)).factor().is_positive
        assert ConicNormFactor(circle) == 1
        assert ConicNormFactor(-circle) == -1

    def test_zero_radius_circle(self):
        center = symbols("x y", real=True)
        circle = Circle(center, 0)
        assert QuadraticForm(circle, ORIGIN).is_nonpositive
        assert ConicNormFactor(circle) == 1
        assert ConicNormFactor(-circle) == -1

    def test_point_conic(self):
        point = symbols("x y z", positive=True)
        conic = PointConic(point)
        assert QuadraticForm(conic, ORIGIN).is_nonpositive
        assert ConicNormFactor(conic) == 1
        assert ConicNormFactor(-conic) == -1

    def test_line_pair(self):
        line1 = Matrix(symbols("a b c", real=True))
        line2 = Matrix(symbols("d e f", real=True))
        line_pair = LinePair(line1, line2)
        assert ConicNormFactor(line_pair) == 1
        assert ConicNormFactor(-line_pair) == 1

    def test_zero_conic_matrix(self):
        conic = Matrix.zeros(3, 3)
        assert ConicNormFactor(conic) == 1

    def test_undecidable(self):
        conic = ConicMatrix(*symbols("a b c d e f", real=True))
        assert isinstance(ConicNormFactor(conic), ConicNormFactor)

    def test_function_properties(self):
        conic = ConicMatrix(*symbols("a b c d e f", real=True))
        factor = ConicNormFactor(conic)
        assert factor.is_nonzero is True
        assert factor.is_integer is True
        assert factor.is_positive is None

    def test_simplify(self):
        conic = ConicMatrix(*symbols("a b c d e f", real=True))
        f = ConicNormFactor(conic)
        assert Abs(f) == 1
        assert f * f == 1
        assert f**3 == f


class TestIsDegenerate:
    def test_numeric(self):
        assert IsDegenerate(Circle((1, 2), 0)) is True
        assert IsDegenerate(UNIT_CIRCLE) is False

    def test_symbolic_conics(self):
        assert IsDegenerate(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None
        assert IsDegenerate(Matrix.diag(symbols("a,c,f", positive=True))) is False

    def test_symbolic_conic_from_focus_and_directrix(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        eccentricity = symbols("e", positive=True)
        conic = ConicFromFocusAndDirectrix(focus, directrix, eccentricity)
        assert IsDegenerate(conic) is False

    def test_symbolic_line_pair(self):
        line_pair = LinePair(Matrix(symbols("a,b,c")), Matrix(symbols("d,e,f")))
        assert IsDegenerate(line_pair) is True

    def test_symbolic_point_conics(self):
        zero_circle = Circle(symbols("x,y"), 0)
        assert IsDegenerate(zero_circle) is True

        point_conic = PointConic(symbols("x,y,z"))
        assert IsDegenerate(point_conic) is True


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

    def test_symbolic_point_conics(self):
        zero_circle = Circle(symbols("x,y"), 0)
        assert IsFiniteConic(zero_circle) is True

        finite_point_conic = PointConic(symbols("x,y", real=True))
        assert IsFiniteConic(finite_point_conic) is True

        ideal_point_conic = PointConic([*symbols("x,y"), 0])
        assert IsFiniteConic(ideal_point_conic) is False

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

    def test_point_conic(self):
        point_conic = PointConic(symbols("x,y,z"))
        assert IsEllipse(point_conic) is False


class TestIsCircle:
    def test_symbolic_circular_conics(self):
        x, y = symbols("x,y", real=True)

        circle = Circle((x, y), symbols("r", positive=True))
        assert IsCircle(circle) is True

        circle_or_point = Circle((x, y), symbols("r", nonnegative=True))
        assert IsCircle(circle_or_point) is None

        complex_circle = Circle((x, y), symbols("r", positive=True) * I)
        assert IsCircle(complex_circle) is False

        double_ideal_line = LinePair(IDEAL_LINE, IDEAL_LINE)
        assert IsCircle(double_ideal_line) is False

    def test_symbolic_ellipse(self):
        x, y = symbols("x,y", real=True)
        r1, r2 = symbols("r1,r2", positive=True)

        assert IsCircle(Ellipse((x, y), r1, r1)) is True
        assert IsCircle(Ellipse((x, y), r1, r2)) is None
        assert IsCircle(Ellipse((x, y), r1, r1 + r2)) is False


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
        r1_dir = symbols("dx,dy", positive=True)
        assert IsComplexEllipse(Ellipse(center, *r)) is False
        assert IsComplexEllipse(Ellipse(center, *r, r1_direction=r1_dir)) is False
        assert IsComplexEllipse(Ellipse(center, *imag_r)) is True
        assert IsComplexEllipse(Ellipse(center, *imag_r, r1_direction=r1_dir)) is True

    def test_symbolic_point_conics(self):
        zero_circle = Circle(symbols("x,y"), 0)
        assert IsComplexEllipse(zero_circle) is False

        point_conic = PointConic(symbols("x,y,z"))
        assert IsComplexEllipse(point_conic) is False


class TestIsParabola:
    def test_numeric(self):
        assert IsParabola(ConicFromPoly(x * x - y))
        assert IsParabola(LinePair(X_AXIS, HorizontalLine(1))) is False
        assert IsParabola(UNIT_CIRCLE) is False

    def test_symbolic(self):
        assert IsParabola(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_focus_and_directrix(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert IsParabola(parabola) is True

    def test_symbolic_point_conics(self):
        zero_circle = Circle(symbols("x,y"), 0)
        assert IsParabola(zero_circle) is False

        point_conic = PointConic(symbols("x,y,z"))
        assert IsParabola(point_conic) is False


class TestIsHyperbola:
    def test_numeric(self):
        assert IsHyperbola(ConicFromPoly(x * y - 1)) is True
        assert IsHyperbola(ConicFromPoly(x * x - y * y - 1)) is True
        assert IsHyperbola(ConicFromPoly(x * x - y * y)) is False
        assert IsHyperbola(UNIT_CIRCLE) is False

    def test_unspecified(self):
        assert IsHyperbola(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_ellipse_or_parabola(self):
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

    def test_symbolic_point_conics(self):
        zero_circle = Circle(symbols("x,y"), 0)
        assert IsHyperbola(zero_circle) is False

        point_conic = PointConic(symbols("x,y,z"))
        assert IsHyperbola(point_conic) is False


class TestIsCircular:
    def test_circle(self):
        assert IsCircular(UNIT_CIRCLE) is True
        assert IsCircular(Circle((1, 2), 3)) is True

    def test_point_conics(self):
        zero_circle = Circle((1, 2), 0)
        assert IsCircular(zero_circle) is True

        finite_point_conic = PointConic([3, 2, 1])
        assert FocalAxisDirection(finite_point_conic).is_zero_matrix is False
        assert IsCircular(finite_point_conic) is False

        ideal_point_conic = PointConic([*symbols("x,y", positive=True), 0])
        assert IsCircular(ideal_point_conic) is False

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

    def test_symbolic_zero_radius_circle(self):
        r = symbols("r")
        assert IsCircular(Circle((0, 0), r)) is True

    def test_symbolic_ellipse(self):
        x, y = symbols("x,y", real=True)
        r1, r2 = symbols("r1,r2", positive=True)
        assert IsCircular(Ellipse((x, y), r1, r1)) is True
        assert IsCircular(Ellipse((x, y), r1, r2)) is None
        assert IsCircular(Ellipse((x, y), r1, r1 + r2)) is False

    def test_symbolic_conics(self):
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

    def test_symbolic_finite_point_conics(self):
        assert IsLinePair(Circle(symbols("x,y"), 0)) is False
        assert IsLinePair(PointConic(symbols("x,y", real=True))) is False

    def test_symbolic_ideal_point_conic(self):
        ideal_point_conic = PointConic([*symbols("x,y", positive=True), 0])
        assert IsLinePair(ideal_point_conic) is False


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


class TestIsPointConic:
    def test_zero_matrix(self):
        assert IsPointConic(Matrix.zeros(3, 3)) is False

    def test_zero_radius_circle(self):
        circle = Circle(symbols("x,y"), 0)
        assert IsPointConic(circle) is True
        assert IsPointConic(-circle) is True

    def test_finite_point_conic(self):
        point = PointConic(symbols("x,y", real=True))
        assert IsPointConic(point) is True

    def test_ideal_point_conic(self):
        point = PointConic([*symbols("x,y", positive=True), 0])
        assert IsPointConic(point) is True

    def test_line_pair(self):
        line1 = Matrix(symbols("a1,b1,c1", positive=True))
        line2 = Matrix(symbols("a2,b2,c2", positive=True))
        line_pair = LinePair(line1, line2)
        assert IsPointConic(line_pair) is False

    def test_circle(self):
        circle = Circle(symbols("x,y"), symbols("r", positive=True))
        assert IsPointConic(circle) is False

    def test_circle_undecidable(self):
        circle = Circle(symbols("x,y"), symbols("r"))
        assert IsPointConic(circle) is None


class TestIsFinitePointConic:
    def test_zero_radius_circle(self):
        assert IsFinitePointConic(Circle(symbols("x,y"), 0)) is True

    def test_circle(self):
        assert IsFinitePointConic(Circle(symbols("x,y"), 1)) is False

    def test_ideal_point_conic(self):
        ideal_point = ConicFromPoly(x * x + 1)
        assert IsFinitePointConic(ideal_point) is False

    def test_line_pair(self):
        assert IsFinitePointConic(LinePair(X_AXIS, Y_AXIS)) is False
        assert IsFinitePointConic(LinePair(X_AXIS, IDEAL_LINE)) is False
        assert IsFinitePointConic(LinePair(IDEAL_LINE, IDEAL_LINE)) is False

    def test_undecidable(self):
        assert IsFinitePointConic(ConicMatrix(*symbols("a,b,c,d,e,f"))) is None
