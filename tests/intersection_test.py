from sympy import I, Matrix, sqrt, symbols

from lib.circle import UNIT_CIRCLE
from lib.degenerate_conic import ZERO_CONIC, line_pair_conic
from lib.incidence import conic_contains_point
from lib.intersection import conic_x_line, line_x_line
from lib.line import (
    IDEAL_LINE,
    X_AXIS,
    Y_AXIS,
    ZERO_LINE,
    horizontal_line,
    vertical_line,
)
from lib.matrix import NonzeroCross, conic_matrix, is_nonzero_multiple
from lib.point import ORIGIN, ZERO_POINT, ideal_point, point_to_xy


class TestLineXLine:
    def test_horizonal_x_vertical(self):
        x, y = symbols("x,y")
        line1 = vertical_line(x)
        line2 = horizontal_line(y)
        assert is_nonzero_multiple(line_x_line(line1, line2), [x, y, 1])

    def test_parallel_lines(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([1, 2, 4])
        assert line_x_line(line1, line2) == Matrix([2, -1, 0])

    def test_coincident_lines(self):
        line = Matrix([1, 2, 3])
        assert line_x_line(line, line) == ZERO_POINT


class TestConicXLine:
    def test_general_case(self):
        conic = conic_matrix(1, 2, 3, 4, 5, 6)
        line = Matrix([1, 2, 3])
        p1, p2 = conic_x_line(conic, line)
        assert conic_contains_point(conic, p1)
        assert conic_contains_point(conic, p2)

    def test_x_y_axes_conic(self):
        conic = line_pair_conic(X_AXIS, Y_AXIS)
        p1, p2 = conic_x_line(conic, IDEAL_LINE)
        assert is_nonzero_multiple(p1, ideal_point(0, 1))
        assert is_nonzero_multiple(p2, ideal_point(1, 0))

    def test_x_axis_plus_ideal_line_conic(self):
        conic = line_pair_conic(X_AXIS, IDEAL_LINE)
        p1, p2 = conic_x_line(conic, Y_AXIS)
        assert is_nonzero_multiple(p1, ORIGIN)
        assert is_nonzero_multiple(p2, ideal_point(0, 1))

    def test_y_axis_plus_ideal_line_conic(self):
        conic = line_pair_conic(Y_AXIS, IDEAL_LINE)
        p1, p2 = conic_x_line(conic, X_AXIS)
        assert is_nonzero_multiple(p1, ideal_point(1, 0))
        assert is_nonzero_multiple(p2, ORIGIN)

    def test_double_x_axis_conic(self):
        conic = line_pair_conic(X_AXIS, X_AXIS)
        assert is_nonzero_multiple(conic_x_line(conic, Y_AXIS)[0], Matrix([0, 0, 1]))
        assert is_nonzero_multiple(
            conic_x_line(conic, IDEAL_LINE)[0],
            Matrix([1, 0, 0]),
        )

    def test_conic_contains_line(self):
        zeros = (ZERO_POINT, ZERO_POINT)
        assert conic_x_line(line_pair_conic(X_AXIS, Y_AXIS), X_AXIS) == zeros
        assert conic_x_line(line_pair_conic(X_AXIS, X_AXIS), X_AXIS) == zeros

    def test_zero_line(self):
        zeros = (ZERO_POINT, ZERO_POINT)
        assert conic_x_line(line_pair_conic(X_AXIS, Y_AXIS), ZERO_LINE) == zeros

    def test_zero_conic(self):
        zeros = (ZERO_POINT, ZERO_POINT)
        assert conic_x_line(ZERO_CONIC, Matrix([1, 2, 3])) == zeros

    def test_double_intersection(self):
        intersections = conic_x_line(UNIT_CIRCLE, horizontal_line(1))
        assert is_nonzero_multiple(intersections[0], Matrix([0, 1, 1]))
        assert is_nonzero_multiple(intersections[1], Matrix([0, 1, 1]))

    def test_complex_intersection(self):
        intersections = conic_x_line(UNIT_CIRCLE, horizontal_line(2))
        intersections = sorted((tuple(point_to_xy(i)) for i in intersections), key=str)
        assert intersections == sorted([(-sqrt(3) * I, 2), (sqrt(3) * I, 2)], key=str)

    def test_symbolic_conic(self):
        conic = conic_matrix(*symbols("a b c d e f"))
        intersections = conic_x_line(conic, IDEAL_LINE)
        assert isinstance(intersections, NonzeroCross)

    def test_symbolic_line(self):
        line = Matrix(symbols("l1 l2 l3"))
        intersections = conic_x_line(UNIT_CIRCLE, line)
        assert isinstance(intersections, NonzeroCross)

    def test_partially_symbolic_conic(self):
        b, c, d, e, f = symbols("b c d e f")
        conic = conic_matrix(1, b, c, d, e, f)
        p1, p2 = conic_x_line(conic, IDEAL_LINE)
        assert conic_contains_point(conic, p1)
        assert conic_contains_point(conic, p2)
