from sympy import I, Matrix, Poly, Rational, sqrt, symbols
from sympy.abc import x, y

from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import (
    AxisDirection,
    ConicFromFocusAndDirectrix,
    ConicFromPoly,
    ConicThroughPoints,
    Eccentricity,
    IdealPoints,
    SplitToLines,
)
from lib.degenerate_conic import LinePair
from lib.ellipse import Ellipse
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, VerticalLine
from lib.matrix import IsNonZeroMultiple, QuadraticForm
from lib.point import IdealPoint, PointToVec3


def AreProjectiveSetsEqual(set1, set2):
    """Compares two sets of projective points or lines for equality."""
    if len(set1) != len(set2):
        return False

    remaining = list(set2)

    for v1 in set1:
        found = False
        for v2 in remaining:
            if IsNonZeroMultiple(v1, v2):
                remaining.remove(v2)
                found = True
                break
        if not found:
            return False
    return True


class TestConicFromPoly:
    def test_expr(self):
        poly = (x + 2) * (3 * y - 4) + x**2
        point = Matrix([x, y, 1])
        assert poly.equals(QuadraticForm(ConicFromPoly(poly), point))

    def test_poly(self):
        assert ConicFromPoly(Poly(1 - x * x - y * y)) == UNIT_CIRCLE

    def test_custom_variables(self):
        cx, cy = symbols("cx cy")
        assert ConicFromPoly(1 - cx**2 - cy**2, x=cx, y=cy) == UNIT_CIRCLE


class TestConicThroughPoints:
    def test_euclidean_points_no_3_collinear(self):
        p1, p2, p3, p4, p5 = (1, 2), (2, 3), (3, 5), (5, 8), (8, 13)
        conic = ConicThroughPoints(p1, p2, p3, p4, p5)
        for p in p1, p2, p3, p4, p5:
            assert QuadraticForm(conic, PointToVec3(p)) == 0

    def test_line_pair(self):
        p1, p2, p3, p4, p5 = (0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)
        conic = ConicThroughPoints(p1, p2, p3, p4, p5)
        assert IsNonZeroMultiple(conic, LinePair(X_AXIS, Y_AXIS))

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


class TestAxisDirection:
    def test_circle(self):
        assert AxisDirection(UNIT_CIRCLE).is_zero_matrix

    def test_ellipse(self):
        ellipse = Ellipse((6, 5), 4, 3, r1_direction=(2, 1))
        assert IsNonZeroMultiple(AxisDirection(ellipse), (2, 1, 0))
        assert IsNonZeroMultiple(AxisDirection(-ellipse), (2, 1, 0))

    def test_hyperbola(self):
        hyperbola = ConicFromPoly(x * y - 1)
        assert IsNonZeroMultiple(AxisDirection(hyperbola), (1, 1, 0))
        assert IsNonZeroMultiple(AxisDirection(-hyperbola), (1, 1, 0))

    def test_degenerate_conic(self):
        conic = LinePair(X_AXIS, Y_AXIS)
        dir1 = AxisDirection(conic)
        dir2 = AxisDirection(-conic)
        assert AreProjectiveSetsEqual([dir1, dir2], [[1, 1, 0], [1, -1, 0]])


class TestSplitToLines:
    def test_horizontal_lines(self):
        lines = [HorizontalLine(1), HorizontalLine(2)]
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_vertical_lines(self):
        lines = [VerticalLine(1), VerticalLine(2)]
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_double_line(self):
        lines = [Matrix([1, 2, 3])] * 2
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_double_ideal_line(self):
        lines = [IDEAL_LINE] * 2
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_general_case(self):
        lines = [Matrix([1, 2, 3]), Matrix([4, 5, 6])]
        conic = LinePair(*lines)
        assert AreProjectiveSetsEqual(lines, SplitToLines(conic))

    def test_point_conic(self):
        conic = Circle((0, 0), 0)
        assert AreProjectiveSetsEqual(
            [Matrix([1, I, 0]), Matrix([1, -I, 0])], SplitToLines(conic)
        )


class TestIdealPoints:
    def test_xy_hyperbola(self):
        hyperbola = ConicFromPoly(x * y)
        ideal_points = IdealPoints(hyperbola)
        ideal_x = IdealPoint(1, 0)
        ideal_y = IdealPoint(0, 1)
        assert AreProjectiveSetsEqual(ideal_points, [ideal_x, ideal_y])

    def test_unit_hyperbola(self):
        hyperbola = ConicFromPoly(x * x - y * y - 1)
        ideal_points = IdealPoints(hyperbola)
        ideal1 = IdealPoint(1, 1)
        ideal2 = IdealPoint(1, -1)
        assert AreProjectiveSetsEqual(ideal_points, [ideal1, ideal2])

    def test_parabola(self):
        parabola = ConicFromPoly(x * x - y)
        ideal_points = IdealPoints(parabola)
        assert IsNonZeroMultiple(ideal_points[0], IdealPoint(0, 1))
        assert IsNonZeroMultiple(ideal_points[1], IdealPoint(0, 1))

    def test_circle(self):
        ideal_points = IdealPoints(UNIT_CIRCLE)
        assert AreProjectiveSetsEqual(
            ideal_points, [IdealPoint(1, I), IdealPoint(1, -I)]
        )
