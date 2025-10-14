from sympy import (
    AppliedPredicate,
    Expr,
    I,
    Matrix,
    Q,
    Rational,
    nan,
    symbols,
)
from sympy.abc import x, y

from lib.central_conic import (
    ConicCenter,
    ConicFromCenterAndPoints,
    ConicFromFociAndRadius,
    LinearEccentricity,
    PrimaryRadius,
    SecondaryRadius,
    ShrinkConicToZero,
)
from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import ConicFromFocusAndDirectrix, ConicFromPoly, IdealPoints
from lib.conic_classification import IsPointConic
from lib.degenerate_conic import LinePair, PointConic
from lib.ellipse import Ellipse, EllipseFromFociAndPoint
from lib.hyperbola import HyperbolaFromFociAndPoint
from lib.line import X_AXIS, HorizontalLine
from lib.matrix import IsNonZeroMultiple
from lib.point import ORIGIN
from lib.sympy_utils import FactorAbs, FactorRadicals
from lib.transform import ScaleXY, TransformConic
from tests.utils import AreProjectiveSetsEqual


class TestConicFromFociAndRadius:
    def test_circle(self):
        center = (1, 2)
        radius = 3
        conic = ConicFromFociAndRadius(center, center, radius)
        assert IsNonZeroMultiple(conic, Circle(center, radius))

    def test_radius_sign_does_not_matter(self):
        conic1 = ConicFromFociAndRadius((1, 2), (3, 4), 5)
        conic2 = ConicFromFociAndRadius((1, 2), (3, 4), -5)
        assert conic1 == conic2


class TestConicFromCenterAndPoints:
    def test_ambiguous_solution(self):
        center = (4, 2)
        p1 = (3, 2)
        p2 = (5, 2)
        p3 = (4, 3)
        conic = ConicFromCenterAndPoints(center, p1, p2, p3)
        assert conic == Matrix.zeros(3, 3)

    def test_circle_centered_at_origin(self):
        p1 = (3, 4)
        p2 = (4, 3)
        p3 = (5, 0)
        conic = ConicFromCenterAndPoints(ORIGIN, p1, p2, p3)
        assert IsNonZeroMultiple(conic, Circle(ORIGIN, 5))

    def test_translated_circle(self):
        center = (1, 3)
        p1 = (4, 7)
        p2 = (5, 6)
        p3 = (6, 3)
        conic = ConicFromCenterAndPoints(center, p1, p2, p3)
        assert IsNonZeroMultiple(conic, Circle(center, 5))

    def test_hyperbola_centered_at_origin(self):
        p1 = (1, 6)
        p2 = (2, 3)
        p3 = (3, 2)
        conic = ConicFromCenterAndPoints(ORIGIN, p1, p2, p3)
        assert IsNonZeroMultiple(conic, ConicFromPoly(x * y - 6))

    def test_parallel_lines(self):
        center = (1, 2)
        p1 = (2, 4)
        p2 = (3, 4)
        p3 = (4, 0)
        conic = ConicFromCenterAndPoints(center, p1, p2, p3)
        assert IsNonZeroMultiple(conic, LinePair(X_AXIS, HorizontalLine(4)))


class TestConicCenter:
    def test_circle(self):
        x, y, r = symbols("x,y,r")
        circle = Circle((x, y), r)
        assert (x, y) == ConicCenter(circle)


class TestSemiAxisLengths:
    def test_circle_radius(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = Circle(center, r)
        assert r == PrimaryRadius(circle)
        assert r == SecondaryRadius(circle)

    def test_ellipse(self):
        ellipse = TransformConic(UNIT_CIRCLE, ScaleXY(2, 3))
        assert PrimaryRadius(ellipse) == 3
        assert PrimaryRadius(-ellipse) == 3
        assert SecondaryRadius(ellipse) == 2
        assert SecondaryRadius(-ellipse) == 2

    def test_complex_ellipse(self):
        ellipse = Ellipse((1, 2), 3 * I, 4 * I)
        assert PrimaryRadius(ellipse) == 3 * I
        assert PrimaryRadius(-ellipse) == 3 * I
        assert SecondaryRadius(ellipse) == 4 * I
        assert SecondaryRadius(-ellipse) == 4 * I

    def test_complex_ellipse_from_foci_and_radius(self):
        ellipse = ConicFromFociAndRadius((-3, 0), (3, 0), 4 * I)
        assert PrimaryRadius(ellipse) == 4 * I
        assert SecondaryRadius(ellipse) == 5 * I

    def test_symbolic_hyperbola(self):
        directrix = Matrix(symbols("a,b,c", positive=True))
        hyperbola = ConicFromFocusAndDirectrix((0, 0), directrix, 2)
        assert FactorRadicals(PrimaryRadius(hyperbola)).is_real is True
        assert FactorRadicals(PrimaryRadius(-hyperbola)).is_real is True
        assert FactorRadicals(SecondaryRadius(hyperbola)).is_real is False
        assert FactorRadicals(SecondaryRadius(-hyperbola)).is_real is False

    def test_parabola(self):
        parabola = ConicFromPoly(x * x - y)
        assert PrimaryRadius(parabola).is_infinite
        assert PrimaryRadius(-parabola).is_infinite
        assert SecondaryRadius(parabola).is_infinite
        assert SecondaryRadius(-parabola).is_infinite

    def test_line_pair(self):
        line1 = Matrix(symbols("a b c", real=True))
        line2 = Matrix(symbols("d e f", real=True))
        line_pair = LinePair(line1, line2)
        assert PrimaryRadius(line_pair) == 0
        assert PrimaryRadius(line_pair) == 0

    def test_finite_point_conic(self):
        zero_circle = Circle(symbols("x y", real=True), 0)
        assert PrimaryRadius(zero_circle) == 0
        assert SecondaryRadius(zero_circle) == 0

        finite_point_conic = PointConic(symbols("x,y", real=True))
        assert PrimaryRadius(finite_point_conic) == 0
        assert SecondaryRadius(finite_point_conic) == 0

    def test_ideal_point_conic(self):
        ideal_point_conic = PointConic([*symbols("x,y"), 0])
        assert PrimaryRadius(ideal_point_conic) == nan
        assert SecondaryRadius(ideal_point_conic) == nan


class TestLinearEccentricity:
    def test_symbolic_circle(self):
        circle = Circle(symbols("x,y", real=True), symbols("r", real=True))
        assert LinearEccentricity(circle) == 0

    def test_symbolic_complex_circle(self):
        circle = Circle(symbols("x,y", real=True), symbols("r", real=True) * I)
        assert LinearEccentricity(circle) == 0

    def test_symbolic_central_conic_standard_form(self):
        fx, r = symbols("fx r", real=True, nonzero=True)
        conic = ConicFromFociAndRadius((-fx, 0), (fx, 0), r)

        def SimplifiedLinEcc(conic: Matrix, assumptions: AppliedPredicate) -> Expr:
            return (
                FactorAbs(LinearEccentricity(conic).factor())
                .refine(assumptions)
                .cancel(r - fx)
            )

        # ellipse
        assert SimplifiedLinEcc(conic, Q.positive(fx) & Q.positive(r - fx)) == fx
        assert SimplifiedLinEcc(-2 * conic, Q.positive(fx) & Q.positive(r - fx)) == fx
        assert SimplifiedLinEcc(conic, Q.negative(fx) & Q.positive(r + fx)) == -fx

        # hyperbola
        assert SimplifiedLinEcc(conic, Q.positive(r) & Q.positive(fx - r)) == fx
        assert SimplifiedLinEcc(-2 * conic, Q.positive(r) & Q.positive(fx - r)) == fx
        assert SimplifiedLinEcc(conic, Q.negative(fx) & Q.positive(r + fx)) == -fx

    def test_numeric_ellipse(self):
        ellipse = ConicFromFociAndRadius((1, 2), (4, 6), 5)
        assert LinearEccentricity(ellipse) == Rational(5, 2)

    def test_numeric_hyperbola(self):
        hyperbola = ConicFromFociAndRadius((1, 2), (4, 6), 2)
        assert LinearEccentricity(hyperbola) == Rational(5, 2)

    def test_symbolic_parabola(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert LinearEccentricity(parabola).is_infinite

    def test_complex_ellipse(self):
        complex_ellipse = Ellipse((1, 2), 4 * I, 5 * I)
        assert LinearEccentricity(complex_ellipse) == 3

    def test_parallel_lines(self):
        c1, c2 = symbols("c1 c2", real=True)
        line1 = Matrix([1, 2, c1])
        line2 = Matrix([1, 2, c2])
        line_pair = LinePair(line1, line2)
        assert LinearEccentricity(line_pair) == nan

    def test_crossing_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        line_pair = LinePair(line1, line2)
        assert LinearEccentricity(line_pair) == 0


class TestShrinkToZero:
    def test_hyperbola(self):
        hyperbola = HyperbolaFromFociAndPoint((1, 2), (3, 4), (0, 0))
        shrunk = ShrinkConicToZero(hyperbola)
        assert shrunk.det() == 0
        assert ConicCenter(hyperbola) == ConicCenter(shrunk)
        assert AreProjectiveSetsEqual(IdealPoints(hyperbola), IdealPoints(shrunk))

    def test_circle(self):
        circle = Circle(symbols("x y"), symbols("r"))
        shrunk = ShrinkConicToZero(circle)
        assert IsPointConic(shrunk)
        assert ConicCenter(shrunk) == (x, y)

    def test_ellipse(self):
        ellipse = EllipseFromFociAndPoint((1, 2), (3, 4), (0, 0))
        shrunk = ShrinkConicToZero(ellipse)
        assert IsPointConic(shrunk)
        assert ConicCenter(ellipse) == ConicCenter(shrunk)
