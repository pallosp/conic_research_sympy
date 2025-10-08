from sympy import I, Matrix, Poly, Rational, sqrt, symbols
from sympy.abc import x, y

from lib.central_conic import ShrinkConicToZero
from lib.circle import UNIT_CIRCLE, Circle
from lib.conic import (
    ConicContainsLine,
    ConicContainsPoint,
    ConicFromFocusAndDirectrix,
    ConicFromPoly,
    ConicNormFactor,
    ConicThroughPoints,
    Eccentricity,
    FocalAxisDirection,
    IdealPoints,
    PolarLine,
    PolePoint,
    ProjectiveConicCenter,
)
from lib.degenerate_conic import LinePair, PointConic
from lib.ellipse import Ellipse
from lib.intersection import ConicXLine
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, HorizontalLine, LineThroughPoint
from lib.matrix import ConicMatrix, IsNonZeroMultiple, QuadraticForm
from lib.point import ORIGIN, IdealPoint, PointToVec3
from lib.transform import TransformConic, Translate
from tests.utils import AreProjectiveSetsEqual


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
            assert ConicContainsPoint(conic, p)

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
    def test_unit_circle(self):
        assert FocalAxisDirection(UNIT_CIRCLE).is_zero_matrix

    def test_symbolic_circle(self):
        assert FocalAxisDirection(Circle(symbols("x y"), symbols("r"))).is_zero_matrix

    def test_ellipse(self):
        ellipse = Ellipse((6, 5), 4, 3, r1_direction=(2, 1))
        assert IsNonZeroMultiple(FocalAxisDirection(ellipse), (2, 1, 0))
        assert IsNonZeroMultiple(FocalAxisDirection(-ellipse), (2, 1, 0))

    def test_point_conic(self):
        ellipse = Ellipse((6, 5), 4, 3, r1_direction=(2, 1))
        point = ShrinkConicToZero(ellipse)
        assert IsNonZeroMultiple(FocalAxisDirection(point), (2, 1, 0))
        assert IsNonZeroMultiple(FocalAxisDirection(-point), (2, 1, 0))

    def test_hyperbola(self):
        hyperbola = ConicFromPoly(x * y - 1)
        assert IsNonZeroMultiple(FocalAxisDirection(hyperbola), (1, 1, 0))
        assert IsNonZeroMultiple(FocalAxisDirection(-hyperbola), (1, 1, 0))

    def test_parabola(self):
        parabola = ConicFromFocusAndDirectrix((1, 2), Matrix([3, 4, 5]), 1)
        assert IsNonZeroMultiple(FocalAxisDirection(parabola), (3, 4, 0))

    def test_degenerate_conic(self):
        conic = LinePair(X_AXIS, Y_AXIS)
        dir1 = FocalAxisDirection(conic)
        dir2 = FocalAxisDirection(-conic)
        assert AreProjectiveSetsEqual([dir1, dir2], [[1, 1, 0], [1, -1, 0]])


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
            ideal_points,
            [IdealPoint(1, I), IdealPoint(1, -I)],
        )


class TestProjectiveConicCenter:
    def test_circle(self):
        center = symbols("x,y")
        circle = Circle(center, symbols("r"))
        assert ProjectiveConicCenter(circle) == Matrix([*center, 1])

    def test_ellipse(self):
        center = symbols("x,y")
        ellipse = Ellipse(
            center,
            symbols("r1", positive=True),
            symbols("r2", positive=True),
            r1_direction=symbols("dx,dy", positive=True),
        )
        computed_center = ProjectiveConicCenter(ellipse).expand()
        assert IsNonZeroMultiple(computed_center, Matrix([*center, 1]))

    def test_parabola(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        center = ProjectiveConicCenter(parabola)
        ideal_point = IdealPoints(parabola)[0]
        assert IsNonZeroMultiple(center, ideal_point)

    def test_parallel_line_pair(self):
        line1 = HorizontalLine(1)
        line2 = HorizontalLine(2)
        line_pair = LinePair(line1, line2)
        assert ProjectiveConicCenter(line_pair).is_zero_matrix

    def test_euclidean_and_ideal_line(self):
        line = Matrix(symbols("a,b,c", positive=True))
        line_pair = LinePair(line, IDEAL_LINE)
        assert ProjectiveConicCenter(line_pair).is_zero_matrix

    def test_ideal_point_conic(self):
        ideal_point = ConicFromPoly(x * x + 1)
        assert ProjectiveConicCenter(ideal_point).is_zero_matrix


class TestPolePolar:
    def test_pole_polar_reciprocity(self):
        conic = ConicMatrix(*symbols("a b c d e f"))
        pole = Matrix(symbols("x y z"))
        polar = PolarLine(conic, pole)
        assert PolePoint(conic, polar).equals(pole * conic.det())

    def test_polar_of_circle_center(self):
        center = (2, 3)
        circle = Circle(center, 4)
        polar_line = PolarLine(circle, center)
        assert IsNonZeroMultiple(polar_line, IDEAL_LINE)

    def test_polar_of_point_on_conic(self):
        hyperbola = ConicFromPoly(x * y - 6)
        point = Matrix([3, 2, 1])
        polar_line = PolarLine(hyperbola, point)
        assert point.dot(polar_line) == 0
        intersections = ConicXLine(hyperbola, polar_line)
        assert intersections[0] == intersections[1]  # tangent line

    def test_pole_of_line_tangent_to_conic(self):
        circle = Circle((0, 0), 5)
        tangent_line = LineThroughPoint((3, 4), direction=(-4, 3))
        pole_point = PolePoint(circle, tangent_line)
        assert ConicContainsPoint(circle, pole_point)

    def test_line_pair_conic(self):
        conic = LinePair(X_AXIS, Y_AXIS)
        # one of the lines
        assert PolePoint(conic, X_AXIS).is_zero_matrix
        # concurrent line
        assert PolePoint(conic, Matrix([1, 1, 0])).is_zero_matrix
        # any other line
        assert IsNonZeroMultiple(PolePoint(conic, HorizontalLine(1)), ORIGIN)
        assert IsNonZeroMultiple(PolePoint(conic, Matrix([1, 2, 3])), ORIGIN)
        assert IsNonZeroMultiple(PolePoint(conic, IDEAL_LINE), ORIGIN)

    def test_point_conic(self):
        conic = Circle((3, 2), 0)
        assert IsNonZeroMultiple(PolePoint(conic, X_AXIS), (3, 2, 1))
        assert IsNonZeroMultiple(PolePoint(conic, IDEAL_LINE), (3, 2, 1))
        assert PolePoint(conic, HorizontalLine(2)).is_zero_matrix


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
