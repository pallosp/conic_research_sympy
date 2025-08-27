from sympy import Matrix, nan, symbols
from sympy.abc import x, y

from lib.central_conic import (
    ConicCenter,
    ConicFromCenterAndPoints,
    SemiMajorAxis,
    SemiMinorAxis,
)
from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import ConicFromFocusAndDirectrix, ConicFromPoly
from lib.degenerate_conic import LinePair, PointConic
from lib.line import X_AXIS, HorizontalLine
from lib.matrix import IsNonZeroMultiple
from lib.point import ORIGIN
from lib.transform import ScaleXY, TransformConic
from tests.util import FactorRadicals


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
        assert r == SemiMajorAxis(circle)
        assert r == SemiMinorAxis(circle)

    def test_ellipse(self):
        ellipse = TransformConic(UNIT_CIRCLE, ScaleXY(2, 3))
        assert SemiMajorAxis(ellipse) == 3
        assert SemiMajorAxis(ellipse * -1) == 3
        assert SemiMinorAxis(ellipse) == 2
        assert SemiMinorAxis(ellipse * -1) == 2

    def test_symbolic_hyperbola(self):
        directrix = Matrix(symbols("a,b,c", positive=True))
        hyperbola = ConicFromFocusAndDirectrix((0, 0), directrix, 2)
        assert FactorRadicals(SemiMajorAxis(hyperbola)).is_real is True
        assert FactorRadicals(SemiMajorAxis(-hyperbola)).is_real is True
        assert FactorRadicals(SemiMinorAxis(hyperbola)).is_real is False
        assert FactorRadicals(SemiMinorAxis(-hyperbola)).is_real is False

    def test_parabola(self):
        parabola = ConicFromPoly(x * x - y)
        assert SemiMajorAxis(parabola).is_infinite
        assert SemiMajorAxis(-parabola).is_infinite
        assert SemiMinorAxis(parabola).is_infinite
        assert SemiMinorAxis(-parabola).is_infinite

    def test_line_pair(self):
        line1 = Matrix(symbols("a b c"))
        line2 = Matrix(symbols("d e f"))
        line_pair = LinePair(line1, line2)
        assert SemiMajorAxis(line_pair) == 0
        assert SemiMajorAxis(line_pair) == 0

    def test_finite_point_conic(self):
        zero_circle = Circle(symbols("x y"), 0)
        assert SemiMajorAxis(zero_circle) == 0
        assert SemiMinorAxis(zero_circle) == 0

        finite_point_conic = PointConic(symbols("x,y"))
        assert SemiMajorAxis(finite_point_conic) == 0
        assert SemiMinorAxis(finite_point_conic) == 0

    def test_ideal_point_conic(self):
        ideal_point_conic = PointConic([*symbols("x,y"), 0])
        assert SemiMajorAxis(ideal_point_conic) == nan
        assert SemiMinorAxis(ideal_point_conic) == nan
