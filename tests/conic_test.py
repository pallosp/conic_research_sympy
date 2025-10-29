import pytest
from sympy import (
    I,
    Matrix,
    Poly,
    Rational,
    exp,
    log,
    pi,
    sqrt,
    symbols,
)
from sympy.abc import x, y

from lib.central_conic import conic_from_foci_and_radius, shrink_conic_to_zero
from lib.circle import UNIT_CIRCLE, circle
from lib.conic import (
    IdealPoints,
    conic_from_focus_and_directrix,
    conic_from_poly,
    conic_through_points,
    eccentricity,
    focal_axis,
    focal_axis_direction,
    polar_line,
    pole_point,
    projective_conic_center,
)
from lib.degenerate_conic import line_pair_conic
from lib.ellipse import ellipse
from lib.incidence import conic_contains_point
from lib.intersection import conic_x_line
from lib.line import (
    IDEAL_LINE,
    X_AXIS,
    Y_AXIS,
    angle_bisector,
    horizontal_line,
    line_through_point,
    perpendicular_line,
)
from lib.matrix import (
    conic_matrix,
    is_nonzero_multiple,
    is_positive_multiple,
    quadratic_form,
)
from lib.point import ORIGIN, ideal_point, ideal_point_on_line
from lib.transform import rotate, transform_conic
from tests.utils import are_projective_sets_equal


class TestConicFromPoly:
    def test_expr(self):
        poly = (x + 2) * (3 * y - 4) + x**2
        point = Matrix([x, y, 1])
        assert poly.equals(quadratic_form(conic_from_poly(poly), point))

    def test_poly(self):
        assert conic_from_poly(Poly(1 - x * x - y * y)) == UNIT_CIRCLE

    def test_custom_variables(self):
        cx, cy = symbols("cx cy")
        assert conic_from_poly(1 - cx**2 - cy**2, x=cx, y=cy) == UNIT_CIRCLE


class TestConicThroughPoints:
    def test_euclidean_points_no_3_collinear(self):
        p1, p2, p3, p4, p5 = (1, 2), (2, 3), (3, 5), (5, 8), (8, 13)
        conic = conic_through_points(p1, p2, p3, p4, p5)
        for p in p1, p2, p3, p4, p5:
            assert conic_contains_point(conic, p)

    def test_line_pair(self):
        p1, p2, p3, p4, p5 = (0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)
        conic = conic_through_points(p1, p2, p3, p4, p5)
        assert is_nonzero_multiple(conic, line_pair_conic(X_AXIS, Y_AXIS))

    def four_collinear_points(self):
        p1, p2, p3, p4, p5 = (0, 0), (1, 0), (2, 0), (3, 0), (0, 1)
        conic = conic_through_points(p1, p2, p3, p4, p5)
        assert conic.is_zero_matrix

    def coincident_points(self):
        p1, p2, p3, p4, p5 = (0, 0), (0, 0), (1, 0), (0, 1), (1, 1)
        conic = conic_through_points(p1, p2, p3, p4, p5)
        assert conic.is_zero_matrix


class TestEccentricity:
    def test_symbolic_conic_from_focus_and_directrix(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a,b,c", positive=True))
        ecc = symbols("e", positive=True)
        conic = conic_from_focus_and_directrix(focus, directrix, ecc)
        assert ecc == eccentricity(conic).simplify()
        assert ecc == eccentricity(conic * -2).simplify()

    def test_symbolic_imaginary_ellipse_from_focus_and_directrix(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a,b,c", positive=True))
        ecc = symbols("e", positive=True) * I
        ellipse = conic_from_focus_and_directrix(focus, directrix, ecc)
        assert ecc == eccentricity(ellipse).simplify()
        assert ecc == eccentricity(ellipse * -2).simplify()

    def test_symbolic_circle(self):
        center = symbols("x y", real=True)
        radius = symbols("r", real=True)
        assert eccentricity(circle(center, radius)) == 0
        assert eccentricity(circle(center, radius) * -2) == 0

    def test_numeric_parabola(self):
        parabola = conic_from_poly(x * x - y)
        assert eccentricity(parabola) == 1
        assert eccentricity(parabola * -2) == 1

    def test_symbolic_parabola(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a,b,c", positive=True))
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        assert eccentricity(parabola) == 1
        assert eccentricity(parabola * -2) == 1

    def test_numeric_rectangular_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 5)
        assert eccentricity(hyperbola) == sqrt(2)
        assert eccentricity(hyperbola * -2) == sqrt(2)

    def test_numeric_degenerate_conic(self):
        conic = line_pair_conic(X_AXIS, Matrix([24, 7, 0]))
        ecc1, ecc2 = eccentricity(conic), eccentricity(-conic)
        assert sorted([ecc1, ecc2]) == [Rational(5, 4), Rational(5, 3)]

    def test_symbolic_ellipse(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a,b,c", positive=True))
        ecc = Rational(1, 2)
        ellipse = conic_from_focus_and_directrix(focus, directrix, ecc)
        assert eccentricity(ellipse).factor() == ecc
        assert eccentricity(-ellipse).factor() == ecc

    def test_point_conic(self):
        focus = ORIGIN
        directrix = Matrix(symbols("a,b,c", positive=True))
        ecc = Rational(1, 2)
        ellipse = conic_from_focus_and_directrix(focus, directrix, ecc)
        point_conic = shrink_conic_to_zero(ellipse)
        assert eccentricity(point_conic).factor() == ecc
        assert eccentricity(-point_conic).factor() == ecc


class TestAxisDirection:
    def test_unit_circle(self):
        assert focal_axis_direction(UNIT_CIRCLE).is_zero_matrix

    def test_symbolic_circle(self):
        assert focal_axis_direction(circle(symbols("x y"), symbols("r"))).is_zero_matrix

    def test_ellipse(self):
        rotated_ellipse = ellipse((6, 5), 4, 3, r1_direction=(2, 1))
        assert is_nonzero_multiple(focal_axis_direction(rotated_ellipse), (2, 1, 0))
        assert is_nonzero_multiple(focal_axis_direction(-rotated_ellipse), (2, 1, 0))

    def test_imaginary_ellipse(self):
        focus1 = (1, 2)
        focus2 = (3, 4)
        imag_ellipse = conic_from_foci_and_radius(focus1, focus2, 5 * I)
        assert is_nonzero_multiple(focal_axis_direction(imag_ellipse), (1, 1, 0))
        assert is_nonzero_multiple(focal_axis_direction(-imag_ellipse), (1, 1, 0))

    def test_point_conic(self):
        r1_direction = (2, 1, 0)
        rotated_ellipse = ellipse((6, 5), 4, 3, r1_direction=r1_direction)
        point = shrink_conic_to_zero(rotated_ellipse)
        assert is_nonzero_multiple(focal_axis_direction(point), r1_direction)
        assert is_nonzero_multiple(focal_axis_direction(-point), r1_direction)

    def test_hyperbola(self):
        hyperbola = conic_from_poly(x * y - 1)
        assert is_nonzero_multiple(focal_axis_direction(hyperbola), (1, 1, 0))
        assert is_nonzero_multiple(focal_axis_direction(-hyperbola), (1, 1, 0))

    def test_parabola(self):
        parabola = conic_from_focus_and_directrix((1, 2), Matrix([3, 4, 5]), 1)
        assert is_nonzero_multiple(focal_axis_direction(parabola), (3, 4, 0))

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
        assert is_nonzero_multiple(focal_axis_direction(line_pair), (0, 1, 0))
        line_pair = line_pair_conic(X_AXIS, -X_AXIS)
        assert is_nonzero_multiple(focal_axis_direction(line_pair), (1, 0, 0))

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


class TestFocalAxis:
    def test_conic_from_focus_and_directrix(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        ecc = symbols("e", positive=True)
        conic = conic_from_focus_and_directrix(focus, directrix, ecc)
        expected = perpendicular_line(directrix, focus)
        axis = focal_axis(conic).applyfunc(
            lambda coord: coord.expand()
            .factor()
            .rewrite(log)
            .factor(deep=True)
            .rewrite(exp)
            .factor(),
        )
        assert is_nonzero_multiple(axis, expected)

    def test_symbolic_circle(self):
        symbolic_circle = circle(symbols("x y"), symbols("r"))
        assert focal_axis(symbolic_circle).is_zero_matrix

    def test_ellipse(self):
        rotated_ellipse = ellipse((6, 5), 4, 3, r1_direction=(2, 1))
        expected = line_through_point((6, 5), direction=(2, 1))
        assert is_nonzero_multiple(focal_axis(rotated_ellipse), expected)

    def test_zero_ellipse(self):
        r1_direction = symbols("x y", nonzero=True)
        rotated_ellipse = ellipse((4, 3), 2, 1, r1_direction=r1_direction)
        ellipse_axis = focal_axis(rotated_ellipse)
        point_conic = shrink_conic_to_zero(rotated_ellipse)
        point_conic_axis = focal_axis(point_conic)
        assert ellipse_axis == point_conic_axis

    def test_line_pair(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        assert is_nonzero_multiple(
            focal_axis(line_pair_conic(line1, line2)),
            angle_bisector(line1, line2),
        )
        assert is_nonzero_multiple(
            focal_axis(line_pair_conic(line1, -line2)),
            angle_bisector(line1, -line2),
        )

    def test_point_conic(self):
        r1_direction = (2, 1, 0)
        rotated_ellipse = ellipse((6, 5), 4, 3, r1_direction=r1_direction)
        point = shrink_conic_to_zero(rotated_ellipse)
        expected = line_through_point((6, 5), direction=r1_direction)
        assert is_nonzero_multiple(focal_axis(point), expected)
        assert is_nonzero_multiple(focal_axis(-point), expected)


class TestIdealPoints:
    def test_xy_hyperbola(self):
        hyperbola = conic_from_poly(x * y)
        ideal_points = IdealPoints(hyperbola)
        ideal_x = ideal_point(1, 0)
        ideal_y = ideal_point(0, 1)
        assert are_projective_sets_equal(ideal_points, [ideal_x, ideal_y])

    def test_unit_hyperbola(self):
        hyperbola = conic_from_poly(x * x - y * y - 1)
        ideal_points = IdealPoints(hyperbola)
        ideal1 = ideal_point(1, 1)
        ideal2 = ideal_point(1, -1)
        assert are_projective_sets_equal(ideal_points, [ideal1, ideal2])

    def test_parabola(self):
        parabola = conic_from_poly(x * x - y)
        ideal_points = IdealPoints(parabola)
        assert is_nonzero_multiple(ideal_points[0], ideal_point(0, 1))
        assert is_nonzero_multiple(ideal_points[1], ideal_point(0, 1))

    def test_circle(self):
        ideal_points = IdealPoints(UNIT_CIRCLE)
        assert are_projective_sets_equal(
            ideal_points,
            [ideal_point(1, I), ideal_point(1, -I)],
        )


class TestProjectiveConicCenter:
    def test_circle(self):
        center = symbols("x,y")
        symbolic_circle = circle(center, symbols("r"))
        assert projective_conic_center(symbolic_circle) == Matrix([*center, 1])

    def test_ellipse(self):
        center = symbols("x,y")
        symbolic_ellipse = ellipse(
            center,
            symbols("r1", positive=True),
            symbols("r2", positive=True),
            r1_direction=symbols("dx,dy", positive=True),
        )
        computed_center = projective_conic_center(symbolic_ellipse).expand()
        assert is_nonzero_multiple(computed_center, Matrix([*center, 1]))

    def test_parabola(self):
        focus = (0, 0)
        directrix = Matrix(symbols("a,b,c", positive=True))
        parabola = conic_from_focus_and_directrix(focus, directrix, 1)
        center = projective_conic_center(parabola)
        ideal_point = IdealPoints(parabola)[0]
        assert is_nonzero_multiple(center, ideal_point)

    def test_parallel_line_pair(self):
        line1 = horizontal_line(1)
        line2 = horizontal_line(2)
        line_pair = line_pair_conic(line1, line2)
        assert projective_conic_center(line_pair).is_zero_matrix

    def test_euclidean_and_ideal_line(self):
        line = Matrix(symbols("a,b,c", positive=True))
        line_pair = line_pair_conic(line, IDEAL_LINE)
        assert projective_conic_center(line_pair).is_zero_matrix

    def test_ideal_point_conic(self):
        ideal_point = conic_from_poly(x * x + 1)
        assert projective_conic_center(ideal_point).is_zero_matrix


class TestPolePolar:
    def test_pole_polar_reciprocity(self):
        conic = conic_matrix(*symbols("a b c d e f"))
        pole = Matrix(symbols("x y z"))
        polar = polar_line(conic, pole)
        assert pole_point(conic, polar).equals(pole * conic.det())

    def test_polar_of_circle_center(self):
        center = (2, 3)
        conic = circle(center, 4)
        polar = polar_line(conic, center)
        assert is_nonzero_multiple(polar, IDEAL_LINE)

    def test_polar_of_point_on_conic(self):
        hyperbola = conic_from_poly(x * y - 6)
        point = Matrix([3, 2, 1])
        polar = polar_line(hyperbola, point)
        assert point.dot(polar) == 0
        intersections = conic_x_line(hyperbola, polar)
        assert intersections[0] == intersections[1]  # tangent line

    def test_pole_of_line_tangent_to_conic(self):
        conic = circle((0, 0), 5)
        tangent_line = line_through_point((3, 4), direction=(-4, 3))
        pole = pole_point(conic, tangent_line)
        assert conic_contains_point(conic, pole)

    def test_line_pair_conic(self):
        conic = line_pair_conic(X_AXIS, Y_AXIS)
        # one of the lines
        assert pole_point(conic, X_AXIS).is_zero_matrix
        # concurrent line
        assert pole_point(conic, Matrix([1, 1, 0])).is_zero_matrix
        # any other line
        assert is_nonzero_multiple(pole_point(conic, horizontal_line(1)), ORIGIN)
        assert is_nonzero_multiple(pole_point(conic, Matrix([1, 2, 3])), ORIGIN)
        assert is_nonzero_multiple(pole_point(conic, IDEAL_LINE), ORIGIN)

    def test_point_conic(self):
        conic = circle((3, 2), 0)
        assert is_nonzero_multiple(pole_point(conic, X_AXIS), (3, 2, 1))
        assert is_nonzero_multiple(pole_point(conic, IDEAL_LINE), (3, 2, 1))
        assert pole_point(conic, horizontal_line(2)).is_zero_matrix
