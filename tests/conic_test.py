from sympy import I, Matrix, Rational, sqrt, symbols
from sympy.abc import x, y

from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import (
    ConicFromFocusAndDirectrix,
    ConicFromPoly,
    ConicThroughPoints,
    Eccentricity,
    IdealPoints,
)
from lib.degenerate_conic import LinePair
from lib.line import X_AXIS, Y_AXIS
from lib.matrix import IsScalarMultiple, QuadraticForm
from lib.point import IdealPoint, PointToVec3


def test_ConicFromPoly():
    poly = (x + 2) * (3 * y - 4) + x**2
    point = Matrix([x, y, 1])
    assert poly.equals(QuadraticForm(ConicFromPoly(poly), point))


class TestConicThroughPoints:
    def test_euclidean_points_no_3_collinear(self):
        p1, p2, p3, p4, p5 = (1, 2), (2, 3), (3, 5), (5, 8), (8, 13)
        conic = ConicThroughPoints(p1, p2, p3, p4, p5)
        for p in p1, p2, p3, p4, p5:
            assert QuadraticForm(conic, PointToVec3(p)) == 0

    def test_line_pair(self):
        p1, p2, p3, p4, p5 = (0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)
        conic = ConicThroughPoints(p1, p2, p3, p4, p5)
        assert IsScalarMultiple(conic, LinePair(X_AXIS, Y_AXIS))

    def four_collinear_points(self):
        p1, p2, p3, p4, p5 = (0, 0), (1, 0), (2, 0), (3, 0), (0, 1)
        conic = ConicThroughPoints(p1, p2, p3, p4, p5)
        assert conic.is_zero_matrix

    def coincident_points(self):
        p1, p2, p3, p4, p5 = (0, 0), (0, 0), (1, 0), (0, 1), (1, 1)
        conic = ConicThroughPoints(p1, p2, p3, p4, p5)
        assert conic.is_zero_matrix


class TestEccentricity:
    def test_symbolic(self):
        a, b, c, fx, fy, e = symbols("a,b,c,fx,fy,e", nonnegative=True)
        conic = ConicFromFocusAndDirectrix((fx, fy), Matrix([a, b, c]), e)
        assert e == Eccentricity(conic).simplify()

    def test_circle(self):
        assert Eccentricity(Circle((1, 2), 3)) == 0
        assert Eccentricity(Circle((1, 2), 3) * -2) == 0

    def test_parabola(self):
        parabola = ConicFromPoly(x * x - y)
        assert Eccentricity(parabola) == 1
        assert Eccentricity(parabola * -2) == 1

    def test_rectangular_hyperbola(self):
        hyperbola = ConicFromPoly(x * y - 5)
        assert Eccentricity(hyperbola) == sqrt(2)
        assert Eccentricity(hyperbola * -2) == sqrt(2)

    def test_degenerate_conic(self):
        conic = LinePair(X_AXIS, Matrix([24, 7, 0]))
        ecc1, ecc2 = Eccentricity(conic), Eccentricity(-conic)
        assert sorted([ecc1, ecc2]) == [Rational(5, 4), Rational(5, 3)]


class TestIdealPoints:
    def test_xy_hyperbola(self):
        hyperbola = ConicFromPoly(x * y)
        ideal_points = IdealPoints(hyperbola)
        ideal_x = IdealPoint(1, 0)
        ideal_y = IdealPoint(0, 1)
        assert (
            IsScalarMultiple(ideal_points[0], ideal_x)
            and IsScalarMultiple(ideal_points[1], ideal_y)
            or IsScalarMultiple(ideal_points[0], ideal_y)
            and IsScalarMultiple(ideal_points[1], ideal_x)
        )

    def test_unit_hyperbola(self):
        hyperbola = ConicFromPoly(x * x - y * y - 1)
        ideal_points = IdealPoints(hyperbola)
        ideal1 = IdealPoint(1, 1)
        ideal2 = IdealPoint(1, -1)
        assert (
            IsScalarMultiple(ideal_points[0], ideal1)
            and IsScalarMultiple(ideal_points[1], ideal2)
            or IsScalarMultiple(ideal_points[0], ideal2)
            and IsScalarMultiple(ideal_points[1], ideal1)
        )

    def test_parabola(self):
        parabola = ConicFromPoly(x * x - y)
        ideal_points = IdealPoints(parabola)
        assert IsScalarMultiple(ideal_points[0], IdealPoint(0, 1))
        assert IsScalarMultiple(ideal_points[1], IdealPoint(0, 1))

    def test_circle(self):
        ideal_points = IdealPoints(UNIT_CIRCLE)
        assert (
            IsScalarMultiple(ideal_points[0], IdealPoint(1, I))
            and IsScalarMultiple(ideal_points[1], IdealPoint(1, -I))
            or IsScalarMultiple(ideal_points[0], IdealPoint(1, -I))
            and IsScalarMultiple(ideal_points[1], IdealPoint(1, I))
        )
