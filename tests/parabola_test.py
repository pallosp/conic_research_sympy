import pytest
from sympy import Matrix, factor, gcd, symbols

from lib.circle import UNIT_CIRCLE
from lib.conic import ConicFromFocusAndDirectrix, PolarLine
from lib.degenerate_conic import LinePair, PointConic
from lib.distance import PointLineDistance
from lib.line import IDEAL_LINE, ArePerpendicular, LineContainsPoint
from lib.matrix import ConicMatrix, IsNonZeroMultiple
from lib.parabola import (
    ParabolaAxis,
    ParabolaDirectrix,
    ParabolaFocalParameter,
    ParabolaFocus,
    ParabolaVertex,
)
from lib.point import ORIGIN, Centroid, PerpendicularFoot


class TestParabolaFocusAndDirectrix:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert IsNonZeroMultiple(ParabolaFocus(parabola), (*focus, 1))
        assert IsNonZeroMultiple(ParabolaDirectrix(parabola), directrix)

    def test_numeric_parallel_line_pair(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        line_pair = LinePair(line1, line2)
        assert ParabolaFocus(line_pair).is_zero_matrix
        assert IsNonZeroMultiple(ParabolaDirectrix(line_pair), IDEAL_LINE)

    def test_numeric_double_line(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 3])
        line_pair = LinePair(line1, line2)
        assert ParabolaFocus(line_pair).is_zero_matrix
        assert ParabolaDirectrix(line_pair).is_zero_matrix

    def test_numeric_double_ideal_line(self):
        line_pair = LinePair(IDEAL_LINE, IDEAL_LINE)
        assert ParabolaFocus(line_pair).is_zero_matrix
        assert ParabolaDirectrix(line_pair).is_zero_matrix

    def test_numeric_real_plus_ideal_line(self):
        line_pair = LinePair(Matrix([1, 2, 3]), IDEAL_LINE)
        assert ParabolaFocus(line_pair).is_zero_matrix
        assert IsNonZeroMultiple(ParabolaDirectrix(line_pair), IDEAL_LINE)

    def test_numeric_ideal_point_conic(self):
        conic = PointConic([1, 2, 0])
        assert ParabolaFocus(conic).is_zero_matrix
        assert IsNonZeroMultiple(ParabolaDirectrix(conic), IDEAL_LINE)

    def test_symbolic_parabola_focus_at_origin(self):
        a, b, c = symbols("a b c", positive=True)
        focus = (0, 0)
        directrix = Matrix([a, b, c])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert IsNonZeroMultiple(ParabolaFocus(parabola), ORIGIN)
        assert IsNonZeroMultiple(ParabolaDirectrix(parabola), directrix)

    def test_not_a_parabola(self):
        focus = (1, 2)
        directrix = Matrix([3, 4, 5])
        eccentricity = 6
        hyperbola = ConicFromFocusAndDirectrix(focus, directrix, eccentricity)
        with pytest.raises(ValueError, match="Not a parabola"):
            ParabolaFocus(hyperbola)
        with pytest.raises(ValueError, match="Not a parabola"):
            ParabolaDirectrix(hyperbola)

    def test_pole_polar_relationship(self):
        conic = ConicMatrix(*symbols("a b c d e f"))
        focus = ParabolaFocus(conic)
        directrix = ParabolaDirectrix(conic)
        polar = PolarLine(conic, focus).applyfunc(factor)
        assert directrix == polar / gcd(list(polar))


class TestParabolaVertex:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        foot = PerpendicularFoot(focus, directrix)
        vertex = Centroid(focus, foot)
        assert vertex == ParabolaVertex(parabola)


class TestParabolaAxis:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        axis = ParabolaAxis(parabola)
        assert LineContainsPoint(axis, ParabolaFocus(parabola))
        assert ArePerpendicular(axis, ParabolaDirectrix(parabola))

    def test_not_a_parabola(self):
        with pytest.raises(ValueError, match="Not a parabola"):
            ParabolaAxis(UNIT_CIRCLE)

    def test_double_finite_line(self):
        line = Matrix(symbols("a b c"))
        line_pair = LinePair(line, line)
        assert ParabolaAxis(line_pair).is_zero_matrix

    def test_parallel_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        line_pair = LinePair(line1, line2)
        assert ParabolaAxis(line_pair).is_zero_matrix

    def test_finite_plus_ideal_line(self):
        finite_line = Matrix(symbols("a b c"))
        line_pair = LinePair(finite_line, IDEAL_LINE)
        assert ParabolaAxis(line_pair).is_zero_matrix


class TestParabolaFocalParameter:
    def test_focal_parameter(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        fp = ParabolaFocalParameter(parabola)
        assert fp == PointLineDistance(focus, directrix)
