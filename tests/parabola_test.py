import pytest
from sympy import Matrix, symbols

from lib.conic import ConicFromFocusAndDirectrix
from lib.degenerate_conic import LinePair
from lib.line import IDEAL_LINE
from lib.matrix import IsNonZeroMultiple
from lib.parabola import ParabolaDirectrix


class TestParabolaDirectrix:

    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert IsNonZeroMultiple(ParabolaDirectrix(parabola), directrix)

    def test_numeric_parallel_line_pair(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        line_pair = LinePair(line1, line2)
        assert IsNonZeroMultiple(ParabolaDirectrix(line_pair), IDEAL_LINE)

    def test_numeric_double_line(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 3])
        line_pair = LinePair(line1, line2)
        assert ParabolaDirectrix(line_pair).is_zero_matrix

    def test_numeric_double_ideal_line(self):
        line_pair = LinePair(IDEAL_LINE, IDEAL_LINE)
        assert ParabolaDirectrix(line_pair).is_zero_matrix

    def test_numeric_real_plus_ideal_line(self):
        line_pair = LinePair(Matrix([1, 2, 3]), IDEAL_LINE)
        assert IsNonZeroMultiple(ParabolaDirectrix(line_pair), IDEAL_LINE)

    def test_symbolic_parabola_focus_at_origin(self):
        a, b, c = symbols("a b c", nonzero=True)
        focus = (0, 0)
        directrix = Matrix([a, b, c])
        parabola = ConicFromFocusAndDirectrix(focus, directrix, 1)
        assert IsNonZeroMultiple(ParabolaDirectrix(parabola), directrix)

    def test_not_a_parabola(self):
        focus = (1, 2)
        directrix = Matrix([3, 4, 5])
        eccentricity = 6
        hyperbola = ConicFromFocusAndDirectrix(focus, directrix, eccentricity)
        with pytest.raises(ValueError, match="Not a parabola"):
            ParabolaDirectrix(hyperbola)
