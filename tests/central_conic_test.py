import pytest
from sympy import AppliedPredicate, Expr, I, Matrix, Q, Rational, nan, pi, symbols, zoo
from sympy.abc import x, y

from lib.central_conic import (
    center_to_focus_vector,
    center_to_vertex_vector,
    conic_center,
    conic_from_center_and_points,
    conic_from_foci_and_radius,
    linear_eccentricity,
    primary_radius,
    radius_in_direction,
    secondary_radius,
    shrink_conic_to_zero,
)
from lib.circle import UNIT_CIRCLE, circle
from lib.conic import IdealPoints, conic_from_focus_and_directrix, conic_from_poly
from lib.conic_classes import is_point_conic
from lib.degenerate_conic import line_pair_conic, point_conic
from lib.ellipse import ellipse, ellipse_from_foci_and_point
from lib.hyperbola import hyperbola_from_foci_and_point
from lib.line import X_AXIS, horizontal_line
from lib.matrix import is_nonzero_multiple
from lib.point import ORIGIN
from lib.transform import scale_xy, transform_conic
from tests.utils import are_projective_sets_equal


class TestConicFromFociAndRadius:
    def test_circle(self):
        center = (1, 2)
        radius = 3
        conic = conic_from_foci_and_radius(center, center, radius)
        assert is_nonzero_multiple(conic, circle(center, radius))

    def test_radius_sign_does_not_matter(self):
        conic1 = conic_from_foci_and_radius((1, 2), (3, 4), 5)
        conic2 = conic_from_foci_and_radius((1, 2), (3, 4), -5)
        assert conic1 == conic2


class TestConicFromCenterAndPoints:
    def test_ambiguous_solution(self):
        center = (4, 2)
        p1 = (3, 2)
        p2 = (5, 2)
        p3 = (4, 3)
        conic = conic_from_center_and_points(center, p1, p2, p3)
        assert conic == Matrix.zeros(3, 3)

    def test_circle_centered_at_origin(self):
        p1 = (3, 4)
        p2 = (4, 3)
        p3 = (5, 0)
        conic = conic_from_center_and_points(ORIGIN, p1, p2, p3)
        assert is_nonzero_multiple(conic, circle(ORIGIN, 5))

    def test_translated_circle(self):
        center = (1, 3)
        p1 = (4, 7)
        p2 = (5, 6)
        p3 = (6, 3)
        conic = conic_from_center_and_points(center, p1, p2, p3)
        assert is_nonzero_multiple(conic, circle(center, 5))

    def test_hyperbola_centered_at_origin(self):
        p1 = (1, 6)
        p2 = (2, 3)
        p3 = (3, 2)
        conic = conic_from_center_and_points(ORIGIN, p1, p2, p3)
        assert is_nonzero_multiple(conic, conic_from_poly(x * y - 6))

    def test_parallel_lines(self):
        center = (1, 2)
        p1 = (2, 4)
        p2 = (3, 4)
        p3 = (4, 0)
        conic = conic_from_center_and_points(center, p1, p2, p3)
        assert is_nonzero_multiple(conic, line_pair_conic(X_AXIS, horizontal_line(4)))


class TestConicCenter:
    def test_circle(self):
        x, y, r = symbols("x,y,r")
        symbolic_circle = circle((x, y), r)
        assert conic_center(symbolic_circle) == Matrix([x, y])


class TestSemiAxisLengths:
    def test_circle_radius(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        symbolic_circle = circle(center, r)
        assert r == primary_radius(symbolic_circle)
        assert r == secondary_radius(symbolic_circle)

    def test_ellipse(self):
        ellipse = transform_conic(UNIT_CIRCLE, scale_xy(2, 3))
        assert primary_radius(ellipse) == 3
        assert primary_radius(-ellipse) == 3
        assert secondary_radius(ellipse) == 2
        assert secondary_radius(-ellipse) == 2

    def test_imaginary_ellipse(self):
        an_ellipse = ellipse((1, 2), 3 * I, 4 * I)
        assert primary_radius(an_ellipse) == 3 * I
        assert primary_radius(-an_ellipse) == 3 * I
        assert secondary_radius(an_ellipse) == 4 * I
        assert secondary_radius(-an_ellipse) == 4 * I

    def test_imaginary_ellipse_from_foci_and_radius(self):
        ellipse = conic_from_foci_and_radius((-3, 0), (3, 0), 4 * I)
        assert primary_radius(ellipse) == 4 * I
        assert secondary_radius(ellipse) == 5 * I

    def test_symbolic_hyperbola(self):
        directrix = Matrix(symbols("a,b,c", positive=True))
        hyperbola = conic_from_focus_and_directrix((0, 0), directrix, 2)
        assert primary_radius(hyperbola).factor(deep=True).is_real is True
        assert primary_radius(-hyperbola).factor(deep=True).is_real is True
        assert secondary_radius(hyperbola).factor(deep=True).is_real is False
        assert secondary_radius(-hyperbola).factor(deep=True).is_real is False

    def test_parabola(self):
        parabola = conic_from_poly(x * x - y)
        assert primary_radius(parabola).is_infinite
        assert primary_radius(-parabola).is_infinite
        assert secondary_radius(parabola).is_infinite
        assert secondary_radius(-parabola).is_infinite

    def test_line_pair(self):
        line1 = Matrix(symbols("a b c", real=True))
        line2 = Matrix(symbols("d e f", real=True))
        line_pair = line_pair_conic(line1, line2)
        assert primary_radius(line_pair) == 0
        assert primary_radius(line_pair) == 0

    def test_finite_point_conic(self):
        zero_circle = circle(symbols("x y", real=True), 0)
        assert primary_radius(zero_circle) == 0
        assert secondary_radius(zero_circle) == 0

        finite_point_conic = point_conic(symbols("x,y", real=True))
        assert primary_radius(finite_point_conic) == 0
        assert secondary_radius(finite_point_conic) == 0

    def test_ideal_point_conic(self):
        ideal_point_conic = point_conic([*symbols("x,y"), 0])
        assert primary_radius(ideal_point_conic) == nan
        assert secondary_radius(ideal_point_conic) == nan


class TestRadiusInDirection:
    def test_bad_arguments(self):
        with pytest.raises(ValueError, match="exactly one of"):
            radius_in_direction(UNIT_CIRCLE)
        with pytest.raises(ValueError, match="exactly one of"):
            radius_in_direction(UNIT_CIRCLE, angle=0, direction=(1, 0))
        with pytest.raises(ValueError, match="Invalid direction"):
            radius_in_direction(UNIT_CIRCLE, direction=(1, 2, 3))

    def test_symbolic_circle(self):
        r = symbols("r", positive=True)
        symbolic_circle = circle(symbols("x,y", real=True), r)
        direction = symbols("dx dy", real=True)
        computed_radius = radius_in_direction(symbolic_circle, direction=direction)
        assert r == computed_radius.simplify()

    def test_vertical_ellipse(self):
        v_ellipse = ellipse((1, 2), 3, 4)
        assert radius_in_direction(v_ellipse, direction=Matrix([1, 0])) == 3
        assert radius_in_direction(v_ellipse, direction=Matrix([0, 1, 0])) == 4

    def test_imaginary_ellipse(self):
        im_ellipse = ellipse((1, 2), 3 * I, 4 * I)
        assert radius_in_direction(im_ellipse, direction=(1, 0)) == 3 * I
        assert radius_in_direction(im_ellipse, direction=(0, 1, 0)) == 4 * I

    def test_radius_at_angle(self):
        vertical_ellipse = ellipse((1, 2), 3, 4)
        assert radius_in_direction(vertical_ellipse, angle=0) == 3
        assert radius_in_direction(vertical_ellipse, angle=pi / 2) == 4


class TestLinearEccentricity:
    def test_symbolic_circle(self):
        symbolic_circle = circle(symbols("x,y", real=True), symbols("r", real=True))
        assert linear_eccentricity(symbolic_circle) == 0

    def test_symbolic_imaginary_circle(self):
        center = symbols("x,y", real=True)
        imaginary_radius = symbols("r", real=True) * I
        imaginary_circle = circle(center, imaginary_radius)
        assert linear_eccentricity(imaginary_circle) == 0

    def test_symbolic_central_conic_standard_form(self):
        fx, r = symbols("fx r", real=True, nonzero=True)
        conic = conic_from_foci_and_radius((-fx, 0), (fx, 0), r)

        def simplified_lin_ecc(conic: Matrix, assumptions: AppliedPredicate) -> Expr:
            return (
                linear_eccentricity(conic)
                .factor(deep=True)
                .refine(assumptions)
                .cancel(r - fx)
            )

        # ellipse
        assert simplified_lin_ecc(conic, Q.positive(fx) & Q.positive(r - fx)) == fx
        assert simplified_lin_ecc(-2 * conic, Q.positive(fx) & Q.positive(r - fx)) == fx
        assert simplified_lin_ecc(conic, Q.negative(fx) & Q.positive(r + fx)) == -fx

        # hyperbola
        assert simplified_lin_ecc(conic, Q.positive(r) & Q.positive(fx - r)) == fx
        assert simplified_lin_ecc(-2 * conic, Q.positive(r) & Q.positive(fx - r)) == fx
        assert simplified_lin_ecc(conic, Q.negative(fx) & Q.positive(r + fx)) == -fx

    def test_numeric_ellipse(self):
        ellipse = conic_from_foci_and_radius((1, 2), (4, 6), 5)
        assert linear_eccentricity(ellipse) == Rational(5, 2)

    def test_numeric_hyperbola(self):
        hyperbola = conic_from_foci_and_radius((1, 2), (4, 6), 2)
        assert linear_eccentricity(hyperbola) == Rational(5, 2)

    def test_symbolic_parabola(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert linear_eccentricity(parabola).is_infinite

    def test_imaginary_ellipse(self):
        imaginary_ellipse = ellipse((1, 2), 4 * I, 5 * I)
        assert linear_eccentricity(imaginary_ellipse) == 3

    def test_parallel_lines(self):
        c1, c2 = symbols("c1 c2", real=True)
        line1 = Matrix([1, 2, c1])
        line2 = Matrix([1, 2, c2])
        line_pair = line_pair_conic(line1, line2)
        assert linear_eccentricity(line_pair) == nan

    def test_crossing_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        line_pair = line_pair_conic(line1, line2)
        assert linear_eccentricity(line_pair) == 0


class TestCenterToFocusVector:
    def test_ellipse(self):
        ellipse = conic_from_foci_and_radius((2, 3), (5, 8), 13)
        expected_vec = Matrix([Rational(3, 2), Rational(5, 2)])
        assert center_to_focus_vector(ellipse) == expected_vec
        assert center_to_focus_vector(-ellipse) == expected_vec

    def test_hyperbola(self):
        hyperbola = conic_from_foci_and_radius((13, 8), (5, 3), 2)
        expected_vec = Matrix([-4, -Rational(5, 2)])
        assert center_to_focus_vector(hyperbola) == expected_vec
        assert center_to_focus_vector(-hyperbola) == expected_vec

    def test_parabola(self):
        parabola = conic_from_poly(x * x - y)
        assert center_to_focus_vector(parabola) == Matrix([nan, zoo])
        assert center_to_focus_vector(-parabola) == Matrix([nan, zoo])

    def test_finite_point_conic(self):
        assert center_to_focus_vector(point_conic([1, 2])) == Matrix([0, 0])

    def test_ideal_point_conic(self):
        assert center_to_focus_vector(point_conic([2, 1, 0])) == Matrix([nan, nan])

    def test_crossing_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        line_pair = line_pair_conic(line1, line2)
        assert center_to_focus_vector(line_pair) == Matrix([0, 0])


class TestCenterToVertexVector:
    def test_ellipse(self):
        ellipse = conic_from_foci_and_radius((0, 0), (3, 4), 10)
        expected_vec = Matrix([6, 8])
        assert center_to_vertex_vector(ellipse) == expected_vec
        assert center_to_vertex_vector(-ellipse) == expected_vec

    def test_hyperbola(self):
        hyperbola = conic_from_foci_and_radius((0, 0), (9, 12), 5)
        expected_vec = Matrix([3, 4])
        assert center_to_vertex_vector(hyperbola) == expected_vec
        assert center_to_vertex_vector(-hyperbola) == expected_vec

    def test_parabola(self):
        parabola = conic_from_poly(x * x - y)
        assert center_to_vertex_vector(parabola) == Matrix([nan, zoo])
        assert center_to_vertex_vector(-parabola) == Matrix([nan, zoo])

    def test_finite_point_conic(self):
        assert center_to_vertex_vector(point_conic([1, 2])) == Matrix([0, 0])

    def test_ideal_point_conic(self):
        assert center_to_vertex_vector(point_conic([2, 1, 0])) == Matrix([nan, nan])

    def test_crossing_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        line_pair = line_pair_conic(line1, line2)
        assert center_to_vertex_vector(line_pair) == Matrix([0, 0])


class TestShrinkToZero:
    def test_hyperbola(self):
        hyperbola = hyperbola_from_foci_and_point((1, 2), (3, 4), (0, 0))
        shrunk = shrink_conic_to_zero(hyperbola)
        assert shrunk.det() == 0
        assert conic_center(hyperbola) == conic_center(shrunk)
        assert are_projective_sets_equal(IdealPoints(hyperbola), IdealPoints(shrunk))

    def test_circle(self):
        symbolic_circle = circle(symbols("x y"), symbols("r"))
        shrunk_circle = shrink_conic_to_zero(symbolic_circle)
        assert is_point_conic(shrunk_circle)
        assert conic_center(shrunk_circle) == Matrix([x, y])

    def test_ellipse(self):
        ellipse = ellipse_from_foci_and_point((1, 2), (3, 4), (0, 0))
        shrunk = shrink_conic_to_zero(ellipse)
        assert is_point_conic(shrunk)
        assert conic_center(ellipse) == conic_center(shrunk)
