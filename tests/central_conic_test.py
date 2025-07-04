from sympy import Matrix, symbols
from sympy.abc import x, y

from lib.central_conic import (
    ConicCenter,
    ConicFromCenterAndPoints,
    SemiMajorAxis,
    SemiMinorAxis,
)
from lib.circle import Circle, UNIT_CIRCLE
from lib.conic import ConicFromPoly
from lib.degenerate_conic import LinePair
from lib.line import X_AXIS, HorizontalLine
from lib.matrix import IsNonZeroMultiple
from lib.point import ORIGIN
from lib.transform import ScaleXY, TransformConic


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


class TestAxes:
    def test_circle_radius(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = Circle(center, r)
        assert r == SemiMajorAxis(circle)
        assert r == SemiMinorAxis(circle)

    def test_ellipse_axes(self):
        ellipse = TransformConic(UNIT_CIRCLE, ScaleXY(2, 3))
        assert SemiMajorAxis(ellipse) == 3
        assert SemiMajorAxis(ellipse * -1) == 3
        assert SemiMinorAxis(ellipse) == 2
        assert SemiMinorAxis(ellipse * -1) == 2
