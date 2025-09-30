import pytest
from sympy import Matrix, symbols
from sympy.abc import x, y

from lib.line import (
    IDEAL_LINE,
    X_AXIS,
    Y_AXIS,
    AngleBisector,
    AreParallel,
    ArePerpendicular,
    HorizontalLine,
    LineBetween,
    LineContainsPoint,
    LineThroughPoint,
    ParallelLine,
    PerpendicularBisector,
    PerpendicularLine,
    VerticalLine,
)
from lib.matrix import IsNonZeroMultiple
from lib.point import ORIGIN, Centroid, IdealPointOnLine


class TestLineContainsPoint:
    def test_numeric(self):
        line = Matrix([1, 1, -3])  # x+y=3
        assert LineContainsPoint(line, (2, 1)) is True
        assert LineContainsPoint(line, (2, 2)) is False

    def test_symbolic(self):
        line = Matrix([1, 1, -3])  # x+y=3
        x = symbols("x")
        assert LineContainsPoint(line, (x, 3 - x)) is True
        assert LineContainsPoint(line, (x, 2 - x)) is False
        assert LineContainsPoint(line, (x, 1)) is None


class TestPerpendicularLine:
    def test_numeric(self):
        assert PerpendicularLine(X_AXIS, (0, 1)) == Y_AXIS


class TestLineThroughPoint:
    def test_missing_direction_and_normal(self):
        point = symbols("x,y")
        with pytest.raises(ValueError, match="Specify exactly one"):
            LineThroughPoint(point)

    def test_overspecified_direction_and_normal(self):
        point = symbols("x,y")
        direction = symbols("dx,dy")
        normal = symbols("nx,ny")
        with pytest.raises(ValueError, match="Specify exactly one"):
            LineThroughPoint(point, direction=direction, normal=normal)

    def test_direction_specified(self):
        point = symbols("x,y")
        direction = symbols("dx,dy")
        line = LineThroughPoint(point, direction=direction)
        assert LineContainsPoint(line, point)
        assert IdealPointOnLine(line) == Matrix((*direction, 0))

    def test_normal_vector_specified(self):
        point = symbols("x,y")
        normal = symbols("nx,ny")
        line = LineThroughPoint(point, normal=normal)
        assert LineContainsPoint(line, point)
        assert line[0:2] == list(normal)


class TestAngleBisector:
    def test_bisector_of_axes(self):
        assert X_AXIS.dot(Matrix([1, 1, 1])) > 0
        assert Y_AXIS.dot(Matrix([1, 1, 1])) < 0

        bisector = AngleBisector(X_AXIS, Y_AXIS)
        assert LineContainsPoint(bisector, ORIGIN)
        assert LineContainsPoint(bisector, (1, -1))

        bisector = AngleBisector(X_AXIS, -Y_AXIS)
        assert LineContainsPoint(bisector, ORIGIN)
        assert LineContainsPoint(bisector, (1, 1))

    def test_parallel_lines_same_direction(self):
        line1 = HorizontalLine(1)
        line2 = HorizontalLine(3)
        center_line = AngleBisector(line1, line2)
        assert IsNonZeroMultiple(center_line, IDEAL_LINE)

    def test_parallel_lines_opposite_direction(self):
        line1 = HorizontalLine(1)
        line2 = -HorizontalLine(3)
        center_line = AngleBisector(line1, line2)
        assert IsNonZeroMultiple(center_line, HorizontalLine(2))

    def test_coincident_lines_same_direction(self):
        line = Matrix([1, 2, 3])
        assert AngleBisector(line, line * 2).is_zero_matrix

    def test_coincident_lines_opposite_direction(self):
        line = Matrix([1, 2, 3])
        assert IsNonZeroMultiple(AngleBisector(line, -line), line)

    def test_two_ideal_lines(self):
        assert AngleBisector(IDEAL_LINE, IDEAL_LINE).is_zero_matrix
        assert AngleBisector(IDEAL_LINE, -IDEAL_LINE).is_zero_matrix

    def test_real_and_ideal_lines(self):
        real_line = Matrix([1, 2, 3])
        assert IsNonZeroMultiple(AngleBisector(real_line, IDEAL_LINE), IDEAL_LINE)


class TestPerpendicularBisector:
    def test_symbolic(self):
        x1, y1, x2, y2 = symbols("x1,y1,x2,y2")
        point1 = (x1, y1)
        point2 = (x2, y2)
        bisector = PerpendicularBisector(point1, point2)
        center = Centroid(point1, point2)
        assert LineContainsPoint(bisector, center)
        assert ArePerpendicular(bisector, LineBetween(point1, point2))


class TestAreParallel:
    def test_numeric(self):
        assert AreParallel(HorizontalLine(1), HorizontalLine(2))
        assert AreParallel(HorizontalLine(1), Y_AXIS) is False
        assert AreParallel(Y_AXIS, IDEAL_LINE)
        assert AreParallel(Y_AXIS, ParallelLine(Y_AXIS, (1, 1)))

    def test_symbolic(self):
        assert AreParallel(Matrix([x, x * 2, 0]), Matrix([x * 2, x * 4, 0]))
        assert AreParallel(Matrix([x, x + x, 0]), Matrix([x, x + 1, 0])) is None


class TestArePerpendicular:
    def test_numeric(self):
        assert ArePerpendicular(HorizontalLine(1), Y_AXIS)
        assert ArePerpendicular(HorizontalLine(1), HorizontalLine(2)) is False
        assert ArePerpendicular(HorizontalLine(1), IDEAL_LINE)

    def test_symbolic(self):
        assert ArePerpendicular(HorizontalLine(y), VerticalLine(x))
