from typing import Any

import pytest
from sympy import I, Matrix, nan, pi, sqrt, symbols, zoo
from sympy.abc import x, y

from lib.circle import (
    IMAGINARY_UNIT_CIRCLE,
    UNIT_CIRCLE,
    circle,
    circle_radius,
    director_circle,
)
from lib.conic import conic_from_focus_and_directrix, conic_from_poly
from lib.degenerate_conic import double_line_conic, line_pair_conic, point_conic
from lib.ellipse import ellipse
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS
from lib.matrix import is_nonzero_multiple, quadratic_form


def test_unit_circle():
    assert is_nonzero_multiple(UNIT_CIRCLE, Matrix.diag([1, 1, -1]))


class TestCircle:
    def test_from_center_and_radius(self):
        circle_matrix = circle((1, 2), r=3)
        expected = Matrix([[-1, 0, 1], [0, -1, 2], [1, 2, 4]])
        assert is_nonzero_multiple(circle_matrix, expected)

    def test_from_center_and_point(self):
        assert is_nonzero_multiple(circle((1, 2), point=(4, 2)), circle((1, 2), r=3))
        assert is_nonzero_multiple(circle((1, 2), point=(1, 2, 1)), circle((1, 2), r=0))

    def test_symbolic_center_and_point(self):
        c = circle((x, y), point=(x + 3, y + 4))
        assert is_nonzero_multiple(c, circle((x, y), r=5))

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"r": 3},
            {"point": (4, 2)},
            {"point": (4, 2, -1)},
            {"point": (1, 2)},
        ],
    )
    def test_value_at_center_is_non_negative(self, kwargs: dict[str, Any]):
        center = Matrix([1, 2, 1])
        assert quadratic_form(circle((1, 2), **kwargs), center) >= 0

    def test_point_and_radius_constructions_agree(self):
        assert circle((1, 2), point=(4, 2)) == circle((1, 2), r=3)

    def test_requires_exactly_one_of_radius_and_point(self):
        with pytest.raises(ValueError, match="Exactly one"):
            circle((1, 2))
        with pytest.raises(ValueError, match="Exactly one"):
            circle((1, 2), r=3, point=(4, 2))

    def test_ideal_center(self):
        assert circle((1, 0, 0), r=3).has(nan, zoo)
        assert circle((1, 0, 0), point=(1, 2)).has(nan, zoo)

    def test_ideal_point(self):
        assert is_nonzero_multiple(
            circle((1, 2), point=(1, 3, 0)), double_line_conic(IDEAL_LINE)
        )
        assert circle((1, 0, 0), point=(1, 3, 0)).has(nan, zoo)


class TestCircleRadius:
    def test_real_circle(self):
        assert circle_radius(circle((1, 2), r=3)) == 3
        assert circle_radius(UNIT_CIRCLE * -2) == 1

    def test_imaginary_circle(self):
        assert circle_radius(IMAGINARY_UNIT_CIRCLE) == I
        assert circle_radius(-IMAGINARY_UNIT_CIRCLE) == I
        assert circle_radius(circle((1, 2), r=3 * I)) == 3 * I

    def test_finite_point_conic(self):
        assert circle_radius(circle((1, 2), r=0)) == 0
        assert circle_radius(point_conic((1, 2))) == 0

    def test_symbolic_circular_conic(self):
        center = symbols("x y")
        r = symbols("r", complex=True)
        # √(r²) = { r if 0 ≤ arg(r) < π, -r otherwise }
        assert circle_radius(circle(center, r=r)) == sqrt(r**2)

    def test_symbolic_finite_point_conic(self):
        center = symbols("x y")
        assert circle_radius(point_conic(center)) == 0


class TestDirectorCircle:
    def test_axis_aligned_ellipse(self):
        result = director_circle(circle((1, 2), r=3))
        assert is_nonzero_multiple(result, circle((1, 2), r=3 * sqrt(2)))
        result = director_circle(ellipse((1, 2), 3, 4))
        assert is_nonzero_multiple(result, circle((1, 2), r=5))

    def test_rotated_ellipse(self):
        rotated_ellipse = ellipse((1, 2), 3, 4, r1_angle=pi / 4)
        result = director_circle(rotated_ellipse)
        assert is_nonzero_multiple(result, circle((1, 2), r=5))

    def test_acute_hyperbola(self):
        hyperbola = conic_from_poly(x**2 / 25 - y**2 / 9 - 1)
        assert is_nonzero_multiple(director_circle(hyperbola), circle((0, 0), r=4))

    def test_rectangular_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 1)
        assert is_nonzero_multiple(director_circle(hyperbola), circle((0, 0), r=0))

    def test_obtuse_hyperbola(self):
        hyperbola = conic_from_poly(x**2 / 9 - y**2 / 25 - 1)
        assert is_nonzero_multiple(director_circle(hyperbola), circle((0, 0), r=4 * I))

    def test_imaginary_ellipse(self):
        imaginary_ellipse = conic_from_poly(x**2 / 9 + y**2 / 16 + 1)
        result = director_circle(imaginary_ellipse)
        assert is_nonzero_multiple(result, circle((0, 0), r=5 * I))
        imaginary_circle = circle((1, 2), r=3 * I)
        result = director_circle(imaginary_circle)
        assert is_nonzero_multiple(result, circle((1, 2), r=3 * sqrt(2) * I))

    def test_parabola(self):
        directrix = Matrix([1, 2, 3])
        parabola = conic_from_focus_and_directrix((4, 5), directrix, 1)
        expected = line_pair_conic(directrix, IDEAL_LINE)
        assert is_nonzero_multiple(director_circle(parabola), expected)

    def test_finite_point_conic(self):
        result = director_circle(point_conic((1, 2)))
        assert is_nonzero_multiple(result, circle((1, 2), r=0))

    def test_ideal_point_conic(self):
        result = director_circle(point_conic((1, 2, 0)))
        assert is_nonzero_multiple(result, double_line_conic(IDEAL_LINE))

    def test_intersecting_line_pair(self):
        perpendicular = line_pair_conic(X_AXIS, Matrix([1, 0, -1]))
        assert is_nonzero_multiple(director_circle(perpendicular), circle((1, 0), r=0))
        oblique = conic_from_poly((y - 2 * x) * (y - x))
        assert is_nonzero_multiple(director_circle(oblique), circle((0, 0), r=0))

    def test_parallel_line_pair(self):
        line_pair = line_pair_conic(Y_AXIS, Matrix([1, 0, -1]))
        result = director_circle(line_pair)
        assert is_nonzero_multiple(result, double_line_conic(IDEAL_LINE))

    def test_line_pair_with_ideal_line(self):
        line_pair = line_pair_conic(X_AXIS, IDEAL_LINE)
        result = director_circle(line_pair)
        assert is_nonzero_multiple(result, double_line_conic(IDEAL_LINE))

    def test_double_line(self):
        assert director_circle(double_line_conic(X_AXIS)).is_zero_matrix
        assert director_circle(double_line_conic(IDEAL_LINE)).is_zero_matrix

    def test_zero_matrix(self):
        assert director_circle(Matrix.zeros(3, 3)).is_zero_matrix
