import pytest
from sympy import I, Matrix, Rational, acos, nan, pi, sqrt, symbols
from sympy.abc import x, y

from lib.circle import circle
from lib.conic import IdealPoints, conic_from_poly, projective_conic_center
from lib.conic_classes import UNIT_HYPERBOLA
from lib.conic_direction import focal_axis_direction
from lib.degenerate_conic import line_pair_conic, point_conic
from lib.ellipse import ellipse
from lib.hyperbola import (
    asymptote_conic,
    asymptote_focal_axis_angle,
    hyperbola_from_foci_and_point,
)
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, horizontal_line, line_through_point
from lib.matrix import is_nonzero_multiple, quadratic_form
from lib.point import ORIGIN
from lib.transform import rotate, transform_line, transform_point
from tests.utils import are_projective_sets_equal


class TestHyperbolaFromFociAndPoint:
    def test_numeric(self):
        # https://www.wolframalpha.com/input?i=x*y=4+focus
        f1 = (sqrt(8), sqrt(8))
        f2 = (-sqrt(8), -sqrt(8))
        p = (4, 1)
        hyperbola = hyperbola_from_foci_and_point(f1, f2, p)
        x, y = symbols("x y")
        assert is_nonzero_multiple(hyperbola, conic_from_poly(x * y - 4))


class TestFocalAxisAsymptoteAngle:
    def test_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 1)
        assert asymptote_focal_axis_angle(hyperbola) == pi / 4
        assert asymptote_focal_axis_angle(-hyperbola) == pi / 4

    def test_parabola(self):
        parabola = conic_from_poly(x * x - y)
        assert asymptote_focal_axis_angle(parabola) == 0

    def test_circle(self):
        a_circle = circle((1, 2), 3)
        assert asymptote_focal_axis_angle(a_circle).is_infinite

    def test_ellipse(self):
        standard_ellipse = ellipse(ORIGIN, 5, 3)
        angle = asymptote_focal_axis_angle(standard_ellipse)
        # As of 2025-11-04 angle.is_imaginary evaluates to None.
        # Tracking issue: https://github.com/sympy/sympy/issues/28541
        assert angle.evalf().is_imaginary
        assert angle == acos(Rational(5, 4))

    def test_complex_ellipse(self):
        complex_ellipse = ellipse(ORIGIN, 5 * I, 3 * I)
        angle = asymptote_focal_axis_angle(complex_ellipse)
        assert angle == acos(3 * I / 4)

    def test_crossing_lines(self):
        line = transform_line(X_AXIS, rotate(pi / 3))
        line_pair = line_pair_conic(X_AXIS, line)
        angle = asymptote_focal_axis_angle(line_pair)
        assert angle == pi / 3
        angle = asymptote_focal_axis_angle(-line_pair)
        assert angle == pi / 6

    def test_parallel_lines(self):
        line_pair = line_pair_conic(horizontal_line(1), horizontal_line(2))
        assert asymptote_focal_axis_angle(line_pair) == pi / 2
        assert asymptote_focal_axis_angle(-line_pair) == 0

    def test_conics_containing_the_ideal_line(self):
        line_pair = line_pair_conic(X_AXIS, IDEAL_LINE)
        assert asymptote_focal_axis_angle(line_pair) == nan
        assert asymptote_focal_axis_angle(-line_pair) == nan

        double_ideal_line = line_pair_conic(IDEAL_LINE, IDEAL_LINE)
        assert asymptote_focal_axis_angle(double_ideal_line) == nan

    def test_finite_point_conic(self):
        zero_circle = circle((1, 2), 0)
        assert asymptote_focal_axis_angle(zero_circle).is_infinite

        finite_point_conic = point_conic((1, 2))
        assert asymptote_focal_axis_angle(finite_point_conic).evalf().is_imaginary

    def test_ideal_point_conic(self):
        ideal_point_conic = point_conic((1, 2, 0))
        assert asymptote_focal_axis_angle(ideal_point_conic) == 0


class TestAsymptoteAngleVsIdealPoints:
    @pytest.mark.parametrize(
        "conic",
        [
            conic_from_poly(x * y - 1),
            ellipse(ORIGIN, 1, 2, r1_direction=(3, 4)),
            ellipse(ORIGIN, I, 2 * I, r1_direction=(3, 4)),
            point_conic((1, 2)),
            line_pair_conic(X_AXIS, Matrix([1, 1, -1])),
        ],
    )
    def test_angle_vs_asymptote_conic_consistency(self, conic: Matrix):
        axis_dir = focal_axis_direction(conic)
        rotation = rotate(asymptote_focal_axis_angle(conic))
        center = projective_conic_center(conic)
        asymptote_dir1 = transform_point(axis_dir, rotation)
        asymptote_dir2 = transform_point(axis_dir, rotation.inv())
        asymptotes = [
            line_through_point(center, direction=asymptote_dir1),
            line_through_point(center, direction=asymptote_dir2),
        ]
        assert is_nonzero_multiple(
            line_pair_conic(*asymptotes),
            asymptote_conic(conic),
        )


class TestAsymptoteConic:
    def test_hyperbola(self):
        hyperbola = conic_from_poly(1 - x * y)
        asymptotes = asymptote_conic(hyperbola)
        assert asymptotes == line_pair_conic(X_AXIS, Y_AXIS)


class TestUnitHyperbola:
    def test_formula(self):
        assert conic_from_poly(x * x - y * y - 1) == UNIT_HYPERBOLA

    def test_ideal_points(self):
        ideal_points = [Matrix([1, 1, 0]), Matrix([1, -1, 0])]
        assert are_projective_sets_equal(IdealPoints(UNIT_HYPERBOLA), ideal_points)

    def test_contains_foci(self):
        focus1 = Matrix([sqrt(2), 0, 1])
        focus2 = Matrix([-sqrt(2), 0, 1])
        assert quadratic_form(UNIT_HYPERBOLA, focus1) > 0
        assert quadratic_form(UNIT_HYPERBOLA, focus2) > 0
