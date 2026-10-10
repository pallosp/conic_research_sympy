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
from lib.conic import conic_from_poly
from lib.degenerate_conic import double_line_conic, point_conic
from lib.ellipse import ellipse
from lib.line import IDEAL_LINE
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
        assert director_circle(circle((1, 2), r=3)) == circle((1, 2), r=3 * sqrt(2))
        assert director_circle(ellipse((1, 2), 3, 4)) == circle((1, 2), r=5)

    def test_rotated_ellipse(self):
        rotated_ellipse = ellipse((1, 2), 3, 4, r1_angle=pi / 4)
        assert director_circle(rotated_ellipse) == circle((1, 2), r=5)

    def test_rectangular_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 1)
        assert director_circle(hyperbola) == circle((0, 0), r=0)
