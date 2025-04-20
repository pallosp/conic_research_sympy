from sympy import Matrix, Rational, symbols

from lib.conic import ConicFromFocusAndDirectrix, ConicThroughPoints, Eccentricity
from lib.degenerate_conic import LinePair
from lib.line import X_AXIS, Y_AXIS
from lib.matrix import IsScalarMultiple, QuadraticForm
from lib.point import PointToVec3


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

    def test_degenerate_conic(self):
        conic = LinePair(X_AXIS, Matrix([24, 7, 0]))
        ecc1, ecc2 = Eccentricity(conic), Eccentricity(-conic)
        assert sorted([ecc1, ecc2]) == [Rational(5, 4), Rational(5, 3)]
