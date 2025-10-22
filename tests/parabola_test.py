import pytest
from sympy import Matrix, factor, gcd, symbols

from lib.circle import UNIT_CIRCLE
from lib.conic import conic_from_focus_and_directrix, polar_line
from lib.degenerate_conic import line_pair_conic, point_conic
from lib.distance import point_line_distance
from lib.incidence import line_contains_point
from lib.line import IDEAL_LINE, are_perpendicular, line_normal
from lib.matrix import conic_matrix, is_nonzero_multiple, is_positive_multiple
from lib.parabola import (
    parabola_axis,
    parabola_direction,
    parabola_directrix,
    parabola_focal_parameter,
    parabola_focus,
    parabola_vertex,
)
from lib.point import ORIGIN, centroid, perpendicular_foot
from lib.transform import scale, transform_conic


class TestParabolaFocusAndDirectrix:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert is_nonzero_multiple(parabola_focus(parabola), (*focus, 1))
        assert is_nonzero_multiple(parabola_directrix(parabola), directrix)

    def test_numeric_parallel_line_pair(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        line_pair = line_pair_conic(line1, line2)
        assert parabola_focus(line_pair).is_zero_matrix
        assert is_nonzero_multiple(parabola_directrix(line_pair), IDEAL_LINE)

    def test_numeric_double_line(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 3])
        line_pair = line_pair_conic(line1, line2)
        assert parabola_focus(line_pair).is_zero_matrix
        assert parabola_directrix(line_pair).is_zero_matrix

    def test_numeric_double_ideal_line(self):
        line_pair = line_pair_conic(IDEAL_LINE, IDEAL_LINE)
        assert parabola_focus(line_pair).is_zero_matrix
        assert parabola_directrix(line_pair).is_zero_matrix

    def test_numeric_real_plus_ideal_line(self):
        line_pair = line_pair_conic(Matrix([1, 2, 3]), IDEAL_LINE)
        assert parabola_focus(line_pair).is_zero_matrix
        assert is_nonzero_multiple(parabola_directrix(line_pair), IDEAL_LINE)

    def test_numeric_ideal_point_conic(self):
        conic = point_conic([1, 2, 0])
        assert parabola_focus(conic).is_zero_matrix
        assert is_nonzero_multiple(parabola_directrix(conic), IDEAL_LINE)

    def test_symbolic_parabola_focus_at_origin(self):
        a, b, c = symbols("a b c", positive=True)
        focus = (0, 0)
        directrix = Matrix([a, b, c])
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert is_nonzero_multiple(parabola_focus(parabola), ORIGIN)
        assert is_nonzero_multiple(parabola_directrix(parabola), directrix)

    def test_not_a_parabola(self):
        focus = (1, 2)
        directrix = Matrix([3, 4, 5])
        eccentricity = 6
        hyperbola = conic_from_focus_and_directrix(focus, directrix, eccentricity)
        with pytest.raises(ValueError, match="Not a parabola"):
            parabola_focus(hyperbola)
        with pytest.raises(ValueError, match="Not a parabola"):
            parabola_directrix(hyperbola)

    def test_pole_polar_relationship(self):
        conic = conic_matrix(*symbols("a b c d e f"))
        focus = parabola_focus(conic)
        directrix = parabola_directrix(conic)
        polar = polar_line(conic, focus).applyfunc(factor)
        assert directrix == polar / gcd(list(polar))


class TestParabolaVertex:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        foot = perpendicular_foot(focus, directrix)
        vertex = centroid(focus, foot)
        assert vertex == parabola_vertex(parabola)


class TestParabolaDirection:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        expected_dir = line_normal(directrix, toward=focus)
        assert is_positive_multiple(parabola_direction(parabola), expected_dir)
        assert is_positive_multiple(parabola_direction(-parabola), expected_dir)

        parabola = transform_conic(parabola, scale(-1))
        expected_dir = -expected_dir
        assert is_positive_multiple(parabola_direction(parabola), expected_dir)
        assert is_positive_multiple(parabola_direction(-parabola), expected_dir)

    def test_not_a_parabola(self):
        with pytest.raises(ValueError, match="Not a parabola"):
            parabola_direction(UNIT_CIRCLE)
        with pytest.raises(ValueError, match="Not a parabola"):
            parabola_direction(point_conic((1, 2, 0)))


class TestParabolaAxis:
    def test_numeric_parabola(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        axis = parabola_axis(parabola)
        assert line_contains_point(axis, parabola_focus(parabola))
        assert are_perpendicular(axis, parabola_directrix(parabola))

    def test_not_a_parabola(self):
        with pytest.raises(ValueError, match="Not a parabola"):
            parabola_axis(UNIT_CIRCLE)

    def test_double_finite_line(self):
        line = Matrix(symbols("a b c"))
        line_pair = line_pair_conic(line, line)
        assert parabola_axis(line_pair).is_zero_matrix

    def test_parallel_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        line_pair = line_pair_conic(line1, line2)
        assert parabola_axis(line_pair).is_zero_matrix

    def test_finite_plus_ideal_line(self):
        finite_line = Matrix(symbols("a b c"))
        line_pair = line_pair_conic(finite_line, IDEAL_LINE)
        assert parabola_axis(line_pair).is_zero_matrix


class TestParabolaFocalParameter:
    def test_focal_parameter(self):
        focus = (6, 5)
        directrix = Matrix([4, 3, 2])
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        fp = parabola_focal_parameter(parabola)
        assert fp == point_line_distance(focus, directrix)
