import pytest
from sympy import Matrix, symbols
from sympy.abc import x, y

from lib.incidence import line_contains_point
from lib.line import (
    IDEAL_LINE,
    X_AXIS,
    Y_AXIS,
    angle_bisector,
    are_parallel,
    are_perpendicular,
    horizontal_line,
    line_between,
    line_normal,
    line_through_point,
    parallel_line,
    perpendicular_bisector,
    perpendicular_line,
    vertical_line,
)
from lib.matrix import is_nonzero_multiple
from lib.point import ORIGIN, centroid, ideal_point_on_line


class TestPerpendicularLine:
    def test_numeric(self):
        assert perpendicular_line(X_AXIS, (0, 1)) == Y_AXIS


class TestLineThroughPoint:
    def test_missing_direction_and_normal(self):
        point = symbols("x,y")
        with pytest.raises(ValueError, match="Specify exactly one"):
            line_through_point(point)

    def test_overspecified_direction_and_normal(self):
        point = symbols("x,y")
        direction = symbols("dx,dy")
        normal = symbols("nx,ny")
        with pytest.raises(ValueError, match="Specify exactly one"):
            line_through_point(point, direction=direction, normal=normal)

    def test_direction_specified(self):
        point = symbols("x,y")
        direction = symbols("dx,dy")
        line = line_through_point(point, direction=direction)
        assert line_contains_point(line, point)
        assert ideal_point_on_line(line) == Matrix((*direction, 0))

    def test_normal_vector_specified(self):
        point = symbols("x,y")
        normal = symbols("nx,ny")
        line = line_through_point(point, normal=normal)
        assert line_contains_point(line, point)
        assert line[0:2] == list(normal)


class TestAngleBisector:
    def test_bisector_of_axes(self):
        assert X_AXIS.dot(Matrix([1, 1, 1])) > 0
        assert Y_AXIS.dot(Matrix([1, 1, 1])) < 0

        bisector = angle_bisector(X_AXIS, Y_AXIS)
        assert line_contains_point(bisector, ORIGIN)
        assert line_contains_point(bisector, (1, -1))

        bisector = angle_bisector(X_AXIS, -Y_AXIS)
        assert line_contains_point(bisector, ORIGIN)
        assert line_contains_point(bisector, (1, 1))

    def test_parallel_lines_same_direction(self):
        line1 = horizontal_line(1)
        line2 = horizontal_line(3)
        center_line = angle_bisector(line1, line2)
        assert is_nonzero_multiple(center_line, IDEAL_LINE)

    def test_parallel_lines_opposite_direction(self):
        line1 = horizontal_line(1)
        line2 = -horizontal_line(3)
        center_line = angle_bisector(line1, line2)
        assert is_nonzero_multiple(center_line, horizontal_line(2))

    def test_coincident_lines_same_direction(self):
        line = Matrix([1, 2, 3])
        assert angle_bisector(line, line * 2).is_zero_matrix

    def test_coincident_lines_opposite_direction(self):
        line = Matrix([1, 2, 3])
        assert is_nonzero_multiple(angle_bisector(line, -line), line)

    def test_two_ideal_lines(self):
        assert angle_bisector(IDEAL_LINE, IDEAL_LINE).is_zero_matrix
        assert angle_bisector(IDEAL_LINE, -IDEAL_LINE).is_zero_matrix

    def test_real_and_ideal_lines(self):
        real_line = Matrix([1, 2, 3])
        assert is_nonzero_multiple(angle_bisector(real_line, IDEAL_LINE), IDEAL_LINE)


class TestPerpendicularBisector:
    def test_symbolic(self):
        x1, y1, x2, y2 = symbols("x1,y1,x2,y2")
        point1 = (x1, y1)
        point2 = (x2, y2)
        bisector = perpendicular_bisector(point1, point2)
        center = centroid(point1, point2)
        assert line_contains_point(bisector, center)
        assert are_perpendicular(bisector, line_between(point1, point2))


class TestLineNormal:
    def test_numeric(self):
        line = Matrix([1, 2, 3])
        assert line_normal(line) == Matrix([1, 2, 0])
        assert line_normal(line, toward=(5, 5)) == Matrix([1, 2, 0])
        assert line_normal(line, toward=(-5, -5)) == Matrix([-1, -2, 0])
        assert line_normal(line, toward=(-1, -1)) == Matrix([0, 0, 0])
        assert line_normal(line, toward=(-5, -5, -1)) == Matrix([-1, -2, 0])
        assert line_normal(line, toward=(-1, 0)) == Matrix([1, 2, 0])
        assert line_normal(line, toward=(-1, 0, 0)) == Matrix([-1, -2, 0])


class TestAreParallel:
    def test_numeric(self):
        assert are_parallel(horizontal_line(1), horizontal_line(2))
        assert are_parallel(horizontal_line(1), Y_AXIS) is False
        assert are_parallel(Y_AXIS, IDEAL_LINE)
        assert are_parallel(Y_AXIS, parallel_line(Y_AXIS, (1, 1)))

    def test_symbolic(self):
        assert are_parallel(Matrix([x, x * 2, 0]), Matrix([x * 2, x * 4, 0]))
        assert are_parallel(Matrix([x, x + x, 0]), Matrix([x, x + 1, 0])) is None


class TestArePerpendicular:
    def test_numeric(self):
        assert are_perpendicular(horizontal_line(1), Y_AXIS)
        assert are_perpendicular(horizontal_line(1), horizontal_line(2)) is False
        assert are_perpendicular(horizontal_line(1), IDEAL_LINE)

    def test_symbolic(self):
        assert are_perpendicular(horizontal_line(y), vertical_line(x))
