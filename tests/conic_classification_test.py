from sympy import I, Matrix, Rational, pi, symbols
from sympy.abc import x, y

from lib.circle import IMAGINARY_UNIT_CIRCLE, UNIT_CIRCLE, circle
from lib.conic import conic_from_focus_and_directrix, conic_from_poly
from lib.conic_classes import (
    is_central_conic,
    is_circle,
    is_circular,
    is_degenerate,
    is_double_line,
    is_ellipse,
    is_finite_conic,
    is_finite_point_conic,
    is_hyperbola,
    is_imaginary_ellipse,
    is_line_pair,
    is_nondegenerate,
    is_parabola,
    is_point_conic,
)
from lib.conic_direction import focal_axis_direction
from lib.degenerate_conic import line_pair_conic, point_conic
from lib.ellipse import ellipse, ellipse_from_foci_and_point
from lib.line import IDEAL_LINE, X_AXIS, Y_AXIS, horizontal_line
from lib.matrix import conic_matrix
from lib.point import ORIGIN


class TestIsDegenerate:
    def test_numeric(self):
        assert is_degenerate(circle((1, 2), 0)) is True
        assert is_degenerate(UNIT_CIRCLE) is False

    def test_general_symbolic_conic(self):
        assert is_degenerate(conic_matrix(*symbols("a,b,c,d,e,f"))) is None
        assert is_nondegenerate(conic_matrix(*symbols("a,b,c,d,e,f"))) is None

        assert is_degenerate(Matrix.diag(symbols("a,c,f", positive=True))) is False
        assert is_nondegenerate(Matrix.diag(symbols("a,c,f", positive=True))) is True

    def test_symbolic_conic_from_focus_and_directrix(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        eccentricity = symbols("e", positive=True)
        conic = conic_from_focus_and_directrix(focus, directrix, eccentricity)
        assert is_degenerate(conic) is False

        # focus lies on directrix
        directrix = Matrix([*symbols("a,b"), 0])
        conic = conic_from_focus_and_directrix(focus, directrix, eccentricity)
        assert is_degenerate(conic) is True

    def test_symbolic_conic_from_center_and_radii(self):
        center = symbols("x y")
        r1, r2 = symbols("r1 r2", nonzero=True)
        r1_direction = symbols("dx dy", positive=True)
        conic = ellipse(center, r1, r2, r1_direction=r1_direction)
        assert is_degenerate(conic) is False

        # zero radius
        conic = ellipse(center, r1, 0, r1_direction=r1_direction)
        assert is_degenerate(conic) is True

    def test_symbolic_line_pair(self):
        line1 = Matrix(symbols("a,b,c"))
        line2 = Matrix(symbols("d,e,f"))
        line_pair = line_pair_conic(line1, line2)
        assert is_degenerate(line_pair) is True

    def test_symbolic_point_conics(self):
        zero_circle = circle(symbols("x,y"), 0)
        assert is_degenerate(zero_circle) is True

        symbolic_point_conic = point_conic(symbols("x,y,z"))
        assert is_degenerate(symbolic_point_conic) is True


class TestIsCentralConic:
    def test_symbolic_ellipse_from_center_and_radii(self):
        center = symbols("x,y")
        symbolic_ellipse = ellipse(
            center,
            *symbols("r1,r2", positive=True),
            r1_direction=symbols("dx dy", positive=True),
        )
        assert is_central_conic(symbolic_ellipse) is True

    def test_symbolic_ellipse_from_foci_and_point(self):
        focus1 = ORIGIN
        focus2 = Matrix(symbols("fx fx", positive=True))
        point = Matrix(symbols("x y", negative=True))
        symbolic_ellipse = ellipse_from_foci_and_point(focus1, focus2, point)
        assert is_central_conic(symbolic_ellipse) is True

    def test_symbolic_parabola(self):
        directrix = Matrix(symbols("a b c", positive=True))
        parabola = conic_from_focus_and_directrix(ORIGIN, directrix, eccentricity=1)
        assert is_central_conic(parabola) is False

    def test_symbolic_hyperbola(self):
        directrix = Matrix(symbols("a b c", positive=True))
        ecc = 1 + symbols("e", positive=True)
        hyperbola = conic_from_focus_and_directrix(ORIGIN, directrix, eccentricity=ecc)
        assert is_central_conic(hyperbola) is True

    def test_symbolic_finite_point_conic(self):
        symbolic_point_conic = point_conic(symbols("x y", real=True))
        assert is_central_conic(symbolic_point_conic) is True

    def test_symbolic_ideal_point_conic(self):
        ideal_point_conic = point_conic([*symbols("x y", real=True), 0])
        assert is_central_conic(ideal_point_conic) is False

    def test_symbolic_crossing_lines(self):
        line = Matrix(symbols("a b c", positive=True))
        line_pair = line_pair_conic(line, X_AXIS)
        assert is_central_conic(line_pair) is True

    def test_symbolic_parallel_lines(self):
        a, b = symbols("a b", positive=True)
        c1, c2 = symbols("c1 c2", real=True)
        line_pair = line_pair_conic(Matrix([a, b, c1]), Matrix([a, b, c2]))
        assert is_central_conic(line_pair) is False

    def test_symbolic_finite_plus_ideal_line(self):
        line = Matrix(symbols("a b c", positive=True))
        line_pair = line_pair_conic(line, IDEAL_LINE)
        assert is_central_conic(line_pair) is False

    def test_symbolic_general_line_pair(self):
        line1 = Matrix(symbols("a b c", positive=True))
        line2 = Matrix(symbols("d e f", positive=True))
        line_pair = line_pair_conic(line1, line2)
        assert is_central_conic(line_pair) is None


class TestIsFiniteConic:
    def test_numeric_point(self):
        assert is_finite_conic(circle((1, 2), 0)) is True

    def test_numeric_circle(self):
        assert is_finite_conic(UNIT_CIRCLE) is True
        assert is_finite_conic(circle((1, 2), 3)) is True

    def test_numeric_parabola(self):
        assert is_finite_conic(conic_from_poly(x * x - y)) is False

    def test_numeric_hyperbola(self):
        assert is_finite_conic(conic_from_poly(x * y - 1)) is False

    def test_numeric_line_pair(self):
        assert is_finite_conic(line_pair_conic(X_AXIS, Y_AXIS)) is False

    def test_symbolic_conic(self):
        assert is_finite_conic(conic_matrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_circle(self):
        center = symbols("x,y")
        radius = symbols("r", positive=True)
        symbolic_circle = circle(center, radius)
        assert is_finite_conic(symbolic_circle) is True

    def test_symbolic_imaginary_ellipse(self):
        imaginary_ellipse = Matrix.diag(symbols("a,c,f", positive=True))
        assert is_finite_conic(imaginary_ellipse) is True

    def test_symbolic_point_conics(self):
        zero_circle = circle(symbols("x,y"), 0)
        assert is_finite_conic(zero_circle) is True

        finite_point_conic = point_conic(symbols("x,y", real=True))
        assert is_finite_conic(finite_point_conic) is True

        ideal_point_conic = point_conic([*symbols("x,y"), 0])
        assert is_finite_conic(ideal_point_conic) is False

    def test_symbolic_line_pair(self):
        line1 = Matrix(symbols("a,b,c", real=True))
        line2 = Matrix(symbols("d,e,f", real=True))
        assert is_finite_conic(line_pair_conic(line1, line2)) is False


class TestIsEllipse:
    def test_numeric_circle(self):
        assert is_ellipse(UNIT_CIRCLE) is True
        assert is_ellipse(-UNIT_CIRCLE) is True
        assert is_ellipse(circle((1, 2), 0)) is False

    def test_numeric_ellipse(self):
        assert is_ellipse(ellipse((1, 2), 3, 4)) is True

    def test_numeric_hyperbola(self):
        assert is_ellipse(conic_from_poly(x * y - 1)) is False

    def test_numeric_point(self):
        assert is_ellipse(circle((1, 2), 0)) is False

    def test_imaginary_circle(self):
        assert is_ellipse(IMAGINARY_UNIT_CIRCLE) is False

    def test_undecidable(self):
        assert is_ellipse(conic_matrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_ellipse(self):
        center = symbols("x,y")
        assert is_ellipse(ellipse(center, *symbols("r1,r2", positive=True))) is True
        assert is_ellipse(ellipse(center, *symbols("r1,r2"))) is None
        assert is_ellipse(IMAGINARY_UNIT_CIRCLE) is False
        assert is_ellipse(-IMAGINARY_UNIT_CIRCLE) is False

    def test_point_conic(self):
        symbolic_point_conic = point_conic(symbols("x,y,z"))
        assert is_ellipse(symbolic_point_conic) is False


class TestIsCircle:
    def test_symbolic_circular_conics(self):
        x, y = symbols("x,y", real=True)

        symbolic_circle = circle((x, y), symbols("r", positive=True))
        assert is_circle(symbolic_circle) is True

        circle_or_point = circle((x, y), symbols("r", nonnegative=True))
        assert is_circle(circle_or_point) is None

        imaginary_circle = circle((x, y), symbols("r", positive=True) * I)
        assert is_circle(imaginary_circle) is False

        double_ideal_line = line_pair_conic(IDEAL_LINE, IDEAL_LINE)
        assert is_circle(double_ideal_line) is False

    def test_symbolic_ellipse(self):
        x, y = symbols("x,y", real=True)
        r1, r2 = symbols("r1,r2", positive=True)

        assert is_circle(ellipse((x, y), r1, r1)) is True
        assert is_circle(ellipse((x, y), r1, r2)) is None
        assert is_circle(ellipse((x, y), r1, r1 + r2)) is False


class TestIsImaginaryEllipse:
    def test_imaginary_unit_circle(self):
        assert is_imaginary_ellipse(IMAGINARY_UNIT_CIRCLE)
        assert is_imaginary_ellipse(-IMAGINARY_UNIT_CIRCLE)

    def test_numeric_ellipse(self):
        assert is_imaginary_ellipse(ellipse((1, 2), 3, 4)) is False
        assert is_imaginary_ellipse(ellipse((1, 2), 3 * I, 4 * I)) is True

    def test_symbolic_ellipse(self):
        center = symbols("x,y")
        r = symbols("r1,r2", positive=True)
        imag_r = [r[0] * I, r[1] * I]
        r1_dir = symbols("dx,dy", positive=True)
        assert is_imaginary_ellipse(ellipse(center, *r)) is False
        assert is_imaginary_ellipse(ellipse(center, *r, r1_direction=r1_dir)) is False
        assert is_imaginary_ellipse(ellipse(center, *imag_r)) is True
        assert (
            is_imaginary_ellipse(ellipse(center, *imag_r, r1_direction=r1_dir)) is True
        )

    def test_symbolic_ellipse_from_focus_and_directrix(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a b c", positive=True))
        eccentricity = symbols("e", positive=True) * I
        conic = conic_from_focus_and_directrix(focus, directrix, eccentricity)
        assert is_imaginary_ellipse(conic)

    def test_symbolic_point_conics(self):
        zero_circle = circle(symbols("x,y"), 0)
        assert is_imaginary_ellipse(zero_circle) is False

        symbolic_point_conic = point_conic(symbols("x,y,z"))
        assert is_imaginary_ellipse(symbolic_point_conic) is False


class TestIsParabola:
    def test_numeric(self):
        assert is_parabola(conic_from_poly(x * x - y))
        assert is_parabola(line_pair_conic(X_AXIS, horizontal_line(1))) is False
        assert is_parabola(UNIT_CIRCLE) is False

    def test_symbolic(self):
        assert is_parabola(conic_matrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_focus_and_directrix(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert is_parabola(parabola) is True

    def test_symbolic_point_conics(self):
        zero_circle = circle(symbols("x,y"), 0)
        assert is_parabola(zero_circle) is False

        symbolic_point_conic = point_conic(symbols("x,y,z"))
        assert is_parabola(symbolic_point_conic) is False


class TestIsHyperbola:
    def test_numeric(self):
        assert is_hyperbola(conic_from_poly(x * y - 1)) is True
        assert is_hyperbola(conic_from_poly(x * x - y * y - 1)) is True
        assert is_hyperbola(conic_from_poly(x * x - y * y)) is False
        assert is_hyperbola(UNIT_CIRCLE) is False

    def test_unspecified(self):
        assert is_hyperbola(conic_matrix(*symbols("a,b,c,d,e,f"))) is None

    def test_symbolic_ellipse_or_parabola(self):
        focus = (0, 0)
        # Real line not going through the focus
        directrix = Matrix(symbols("a,b,c", positive=True))

        circle = conic_from_focus_and_directrix(focus, directrix, 0)
        assert is_hyperbola(circle) is False

        ellipse = conic_from_focus_and_directrix(focus, directrix, Rational(1, 2))
        assert is_hyperbola(ellipse) is False

        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert is_hyperbola(parabola) is False

    def test_symbolic_hyperbola(self):
        focus = (0, 0)
        # Real line not going through the focus
        directrix = Matrix(symbols("a,b,c", positive=True))
        # Eccentricity greater than 1
        eccentricity = symbols("e", positive=True) + 1
        hyperbola = conic_from_focus_and_directrix(focus, directrix, eccentricity)
        assert is_hyperbola(hyperbola) is True

    def test_symbolic_point_conics(self):
        zero_circle = circle(symbols("x,y"), 0)
        assert is_hyperbola(zero_circle) is False

        symbolic_point_conic = point_conic(symbols("x,y,z"))
        assert is_hyperbola(symbolic_point_conic) is False


class TestIsCircular:
    def test_circle(self):
        assert is_circular(UNIT_CIRCLE) is True
        assert is_circular(circle((1, 2), 3)) is True

    def test_point_conics(self):
        zero_circle = circle((1, 2), 0)
        assert is_circular(zero_circle) is True

        finite_point_conic = point_conic([3, 2, 1])
        assert focal_axis_direction(finite_point_conic).is_zero_matrix is False
        assert is_circular(finite_point_conic) is False

        ideal_point_conic = point_conic([*symbols("x,y", positive=True), 0])
        assert is_circular(ideal_point_conic) is False

    def test_imaginary_circle(self):
        assert is_circular(circle((1, 2), I)) is True

    def test_ellipse(self):
        assert is_circular(ellipse((0, 0), 1, 1)) is True
        assert is_circular(ellipse((0, 0), 1, 1, r1_angle=pi / 4)) is True
        assert is_circular(ellipse((0, 0), 1, 2)) is False

    def test_hyperbola(self):
        assert is_circular(conic_from_poly(x * y - 1)) is False

    def test_lines(self):
        assert is_circular(line_pair_conic(IDEAL_LINE, IDEAL_LINE)) is False
        assert is_circular(line_pair_conic(IDEAL_LINE, X_AXIS)) is False
        assert is_circular(line_pair_conic(X_AXIS, X_AXIS)) is False
        assert is_circular(line_pair_conic(X_AXIS, Y_AXIS)) is False

    def test_symbolic_zero_radius_circle(self):
        r = symbols("r")
        assert is_circular(circle((0, 0), r)) is True

    def test_symbolic_ellipse(self):
        x, y = symbols("x,y", real=True)
        r1, r2 = symbols("r1,r2", positive=True)
        assert is_circular(ellipse((x, y), r1, r1)) is True
        assert is_circular(ellipse((x, y), r1, r2)) is None
        assert is_circular(ellipse((x, y), r1, r1 + r2)) is False

    def test_symbolic_conics(self):
        pos1, pos2 = symbols("pos1,pos2", positive=True)
        neg = symbols("neg", negative=True)
        assert is_circular(conic_from_poly(pos1 * x**2 + pos1 * y**2)) is True
        assert is_circular(conic_from_poly(pos1 * x**2 + neg * y**2)) is False
        assert is_circular(conic_from_poly(pos1 * x**2 + pos2 * y**2)) is None


class TestIsLinePair:
    def test_zero_matrix(self):
        assert is_line_pair(Matrix.zeros(3, 3)) is False

    def test_numeric_line_pair(self):
        assert is_line_pair(line_pair_conic(X_AXIS, X_AXIS)) is True
        assert is_line_pair(line_pair_conic(X_AXIS, Y_AXIS)) is True
        assert is_line_pair(line_pair_conic(X_AXIS, horizontal_line(1))) is True
        assert is_line_pair(line_pair_conic(X_AXIS, IDEAL_LINE)) is True
        assert is_line_pair(line_pair_conic(IDEAL_LINE, IDEAL_LINE)) is True
        assert is_line_pair(conic_from_poly(x * x - y * y)) is True

    def test_numeric_hyperbola(self):
        assert is_line_pair(conic_from_poly(x * y - 1)) is False

    def test_numeric_parabola(self):
        assert is_line_pair(conic_from_poly(x * x - y)) is False

    def test_numeric_circle(self):
        assert is_line_pair(UNIT_CIRCLE) is False

    def test_undecidable(self):
        assert (
            is_line_pair(conic_matrix(*symbols("a,b,c,d,e,f", positive=True))) is None
        )

    def test_symbolic_line_pair(self):
        line_or_zero_vector = Matrix(symbols("a,b,c", real=True))
        line1 = Matrix(symbols("d,e,f", positive=True))
        line2 = Matrix(symbols("g,h,i", positive=True))
        assert is_line_pair(line_pair_conic(line_or_zero_vector, line1)) is None
        assert is_line_pair(line_pair_conic(line1, line1)) is True
        assert is_line_pair(line_pair_conic(line1, line2)) is True

    def test_symbolic_finite_point_conics(self):
        assert is_line_pair(circle(symbols("x,y"), 0)) is False
        assert is_line_pair(point_conic(symbols("x,y", real=True))) is False

    def test_symbolic_ideal_point_conic(self):
        ideal_point_conic = point_conic([*symbols("x,y", positive=True), 0])
        assert is_line_pair(ideal_point_conic) is False


class TestIsDoubleLine:
    def test_zero_matrix(self):
        assert is_double_line(Matrix.zeros(3, 3)) is False

    def test_numeric_line_pair(self):
        assert is_double_line(line_pair_conic(X_AXIS, X_AXIS)) is True
        assert is_double_line(line_pair_conic(X_AXIS, Y_AXIS)) is False
        assert is_double_line(line_pair_conic(X_AXIS, horizontal_line(1))) is False
        assert is_double_line(line_pair_conic(X_AXIS, IDEAL_LINE)) is False
        assert is_double_line(line_pair_conic(IDEAL_LINE, IDEAL_LINE)) is True

    def test_symbolic_line_pair(self):
        line_or_zeros = Matrix(symbols("a,b,c", real=True))
        real_line1 = Matrix(symbols("d,e,f", positive=True))
        real_line2 = Matrix(symbols("g,h,i", positive=True))
        assert is_double_line(line_pair_conic(real_line1, real_line1)) is True
        assert is_double_line(line_pair_conic(line_or_zeros, line_or_zeros)) is None
        assert is_double_line(line_pair_conic(real_line1, real_line2)) is None
        assert is_double_line(line_pair_conic(real_line1, IDEAL_LINE)) is False

    def test_symbolic_point(self):
        assert is_double_line(circle(symbols("x,y", real=True), 0)) is False


class TestIsPointConic:
    def test_zero_matrix(self):
        assert is_point_conic(Matrix.zeros(3, 3)) is False

    def test_zero_radius_circle(self):
        zero_circle = circle(symbols("x,y"), 0)
        assert is_point_conic(zero_circle) is True
        assert is_point_conic(-zero_circle) is True

    def test_finite_point_conic(self):
        point = point_conic(symbols("x,y", real=True))
        assert is_point_conic(point) is True

    def test_ideal_point_conic(self):
        point = point_conic([*symbols("x,y", positive=True), 0])
        assert is_point_conic(point) is True

    def test_line_pair(self):
        line1 = Matrix(symbols("a b c", positive=True))
        line2 = Matrix(symbols("d e f", positive=True))
        line_pair = line_pair_conic(line1, line2)
        assert is_point_conic(line_pair) is False

    def test_circle(self):
        center = symbols("x,y")
        radius = symbols("r", positive=True)
        symbolic_circle = circle(center, radius)
        assert is_point_conic(symbolic_circle) is False

    def test_circle_undecidable(self):
        center = symbols("x,y")
        radius = symbols("r", nonnegative=True)
        symbolic_circle = circle(center, radius)
        assert is_point_conic(symbolic_circle) is None


class TestIsFinitePointConic:
    def test_zero_radius_circle(self):
        assert is_finite_point_conic(circle(symbols("x,y"), 0)) is True

    def test_circle(self):
        assert is_finite_point_conic(circle(symbols("x,y"), 1)) is False

    def test_ideal_point_conic(self):
        ideal_point = conic_from_poly(x * x + 1)
        assert is_finite_point_conic(ideal_point) is False

    def test_line_pair(self):
        assert is_finite_point_conic(line_pair_conic(X_AXIS, Y_AXIS)) is False
        assert is_finite_point_conic(line_pair_conic(X_AXIS, IDEAL_LINE)) is False
        assert is_finite_point_conic(line_pair_conic(IDEAL_LINE, IDEAL_LINE)) is False

    def test_undecidable(self):
        assert is_finite_point_conic(conic_matrix(*symbols("a,b,c,d,e,f"))) is None
