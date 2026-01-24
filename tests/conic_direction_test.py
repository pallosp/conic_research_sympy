import pytest
from sympy import Abs, I, Matrix, pi, sign, symbols
from sympy.abc import x, y

from lib.central_conic import conic_from_foci_and_radius, shrink_conic_to_zero
from lib.circle import UNIT_CIRCLE, circle
from lib.conic import conic_from_focus_and_directrix, conic_from_poly
from lib.conic_direction import (
    ConicNormFactor,
    conjugate_axis_direction,
    focal_axis_direction,
)
from lib.degenerate_conic import line_pair_conic, point_conic
from lib.ellipse import ellipse
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, angle_bisector
from lib.matrix import (
    conic_matrix,
    is_nonzero_multiple,
    is_positive_multiple,
    quadratic_form,
)
from lib.point import ORIGIN, ideal_point_on_line
from lib.transform import rotate, transform_conic, translate
from tests.utils import are_projective_sets_equal


class TestConicNormFactor:
    def test_parabola_with_focus_at_origin(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert ConicNormFactor(parabola) == 1
        assert ConicNormFactor(-parabola) == -1

    def test_general_parabola(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        translation = translate(symbols("dx dy", real=True))
        parabola = transform_conic(parabola, translation)
        assert ConicNormFactor(parabola) == 1

    def test_symbolic_circle(self):
        center = [*symbols("x y"), 1]
        symbolic_circle = circle(center, symbols("r", positive=True))
        assert quadratic_form(symbolic_circle, Matrix(center)).factor().is_positive
        assert ConicNormFactor(symbolic_circle) == 1
        assert ConicNormFactor(-symbolic_circle) == -1

    def test_zero_radius_circle(self):
        center = symbols("x y", real=True)
        zero_circle = circle(center, 0)
        assert quadratic_form(zero_circle, ORIGIN).is_nonpositive
        assert ConicNormFactor(zero_circle) == 1
        assert ConicNormFactor(-zero_circle) == -1

    def test_point_conic(self):
        point = symbols("x y z", positive=True)
        conic = point_conic(point)
        assert quadratic_form(conic, ORIGIN).is_nonpositive
        assert ConicNormFactor(conic) == 1
        assert ConicNormFactor(-conic) == -1

    def test_line_pair(self):
        line1 = Matrix(symbols("a b c", real=True))
        line2 = Matrix(symbols("d e f", real=True))
        line_pair = line_pair_conic(line1, line2)
        assert ConicNormFactor(line_pair) == 1
        assert ConicNormFactor(-line_pair) == 1

    def test_zero_conic_matrix(self):
        conic = Matrix.zeros(3, 3)
        assert ConicNormFactor(conic) == 1

    def test_undecidable(self):
        conic = conic_matrix(*symbols("a b c d e f", real=True))
        assert isinstance(ConicNormFactor(conic), ConicNormFactor)

    def test_function_properties(self):
        conic = conic_matrix(*symbols("a b c d e f", real=True))
        factor = ConicNormFactor(conic)
        assert factor.is_nonzero is True
        assert factor.is_integer is True
        assert factor.is_positive is None

    def test_simplify_abs(self):
        conic = conic_matrix(*symbols("a b c d e f", real=True))
        assert Abs(ConicNormFactor(conic)) == 1

    def test_simplify_pow(self):
        conic = conic_matrix(*symbols("a b c d e f", real=True))
        f = ConicNormFactor(conic)
        assert f * f == 1
        assert f**3 == f
        assert 1 / f == f

    def test_simplify_nondegenerate_conic(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a,b,c", positive=True))
        eccentricity = symbols("e", positive=True)
        plus_minus_one = (-1) ** symbols("n", integer=True)
        conic = conic_from_focus_and_directrix(focus, directrix, eccentricity)
        conic *= plus_minus_one
        # Ideally, ConicNormFactor(conic) should simplify to (-1)ⁿ,
        # but Sympy 1.14.0 doesn’t recognize that sign((-1)ⁿ) == (-1)ⁿ.
        f = ConicNormFactor(conic)
        assert f == sign(plus_minus_one)
        assert Abs(f) == 1


class TestAxisDirection:
    def test_unit_circle(self):
        assert focal_axis_direction(UNIT_CIRCLE).is_zero_matrix
        assert conjugate_axis_direction(UNIT_CIRCLE).is_zero_matrix

    def test_symbolic_circle(self):
        assert focal_axis_direction(circle(symbols("x y"), symbols("r"))).is_zero_matrix

    def test_ellipse(self):
        h_ellipse = ellipse((5, 4), 3, 2, r1_direction=(1, 0))
        assert is_positive_multiple(focal_axis_direction(h_ellipse), (1, 0, 0))
        assert is_positive_multiple(focal_axis_direction(-h_ellipse), (1, 0, 0))
        assert is_positive_multiple(conjugate_axis_direction(h_ellipse), (0, 1, 0))

        v_ellipse = ellipse((5, 4), 3, 2, r1_direction=(0, 1))
        assert is_positive_multiple(focal_axis_direction(v_ellipse), (0, 1, 0))
        assert is_positive_multiple(focal_axis_direction(-v_ellipse), (0, 1, 0))
        assert is_positive_multiple(conjugate_axis_direction(v_ellipse), (-1, 0, 0))

        rot_ellipse = ellipse((6, 5), 4, 3, r1_direction=(2, 1))
        assert is_positive_multiple(focal_axis_direction(rot_ellipse), (2, 1, 0))
        assert is_positive_multiple(focal_axis_direction(-rot_ellipse), (2, 1, 0))
        assert is_positive_multiple(conjugate_axis_direction(rot_ellipse), (-1, 2, 0))

    def test_imaginary_ellipse(self):
        focus1 = (1, 2)
        focus2 = (3, 4)
        imag_ellipse = conic_from_foci_and_radius(focus1, focus2, 5 * I)
        assert is_positive_multiple(focal_axis_direction(imag_ellipse), (1, 1, 0))
        assert is_positive_multiple(focal_axis_direction(-imag_ellipse), (1, 1, 0))

    def test_point_conic(self):
        r1_direction = (2, 1, 0)
        rotated_ellipse = ellipse((6, 5), 4, 3, r1_direction=r1_direction)
        point = shrink_conic_to_zero(rotated_ellipse)
        assert is_positive_multiple(focal_axis_direction(point), r1_direction)
        assert is_positive_multiple(focal_axis_direction(-point), r1_direction)

    def test_zero_radius_circle(self):
        point = circle((1, 2), 0)
        assert focal_axis_direction(point) == Matrix([0, 0, 0])
        assert conjugate_axis_direction(point) == Matrix([0, 0, 0])

    def test_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 1)
        assert is_positive_multiple(focal_axis_direction(hyperbola), (1, 1, 0))
        assert is_positive_multiple(focal_axis_direction(-hyperbola), (1, 1, 0))

    def test_parabola(self):
        parabola = conic_from_focus_and_directrix((1, 2), Matrix([3, 4, 5]), 1)
        assert is_positive_multiple(focal_axis_direction(parabola), (3, 4, 0))

    def test_angle_range(self):
        right_parabola = conic_from_poly(y**2 - x)
        up_parabola = conic_from_poly(x**2 - y)
        left_parabola = conic_from_poly(y**2 + x)
        down_parabola = conic_from_poly(x**2 + y)
        nw_parabola = transform_conic(up_parabola, rotate(pi / 4))
        sw_parabola = transform_conic(left_parabola, rotate(pi / 4))
        assert is_positive_multiple(focal_axis_direction(right_parabola), (1, 0, 0))
        assert is_positive_multiple(focal_axis_direction(up_parabola), (0, 1, 0))
        assert is_positive_multiple(focal_axis_direction(left_parabola), (1, 0, 0))
        assert is_positive_multiple(focal_axis_direction(down_parabola), (0, 1, 0))
        assert is_positive_multiple(focal_axis_direction(nw_parabola), (1, -1, 0))
        assert is_positive_multiple(focal_axis_direction(sw_parabola), (1, 1, 0))

    def test_line_pair(self):
        conic = line_pair_conic(X_AXIS, Y_AXIS)
        dir1 = focal_axis_direction(conic)
        dir2 = focal_axis_direction(-conic)
        assert are_projective_sets_equal([dir1, dir2], [[1, 1, 0], [1, -1, 0]])

    def test_double_line(self):
        line_pair = line_pair_conic(X_AXIS, X_AXIS)
        assert is_positive_multiple(focal_axis_direction(line_pair), (0, 1, 0))
        line_pair = line_pair_conic(X_AXIS, -X_AXIS)
        assert is_positive_multiple(focal_axis_direction(line_pair), (1, 0, 0))

    @pytest.mark.parametrize(
        ("line1", "line2"),
        [
            (X_AXIS, Y_AXIS),
            (X_AXIS, -Y_AXIS),
            (X_AXIS, IDEAL_LINE),
            (Matrix([1, 2, 3]), Matrix([4, 5, 6])),
        ],
    )
    def test_line_pair_focal_axis_vs_angle_bisector(self, line1: Matrix, line2: Matrix):
        line_pair = line_pair_conic(line1, line2)
        axis_dir = focal_axis_direction(line_pair)
        bisector = angle_bisector(line1, line2)
        assert is_nonzero_multiple(axis_dir, ideal_point_on_line(bisector))
