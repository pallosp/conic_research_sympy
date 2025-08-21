import pytest
from sympy import Matrix, symbols
from sympy.abc import x, y

from lib.line import (
    AreParallel,
    ArePerpendicular,
    HorizontalLine,
    IDEAL_LINE,
    LineBetween,
    LineThroughPoint,
    ParallelLine,
    PerpendicularBisector,
    PerpendicularLine,
    VerticalLine,
    X_AXIS,
    Y_AXIS,
)
from lib.point import IdealPointOnLine, PointToVec3


class TestPerpendicularLine:
    def test_numeric(self):
        assert Y_AXIS == PerpendicularLine(X_AXIS, (0, 1))


class TestLineThroughPoint:
    def test_missing_direction_and_normal(self):
        point = symbols("x,y")
        with pytest.raises(AssertionError):
            LineThroughPoint(point)

    def test_overspecified_direction_and_normal(self):
        point = symbols("x,y")
        direction = symbols("dx,dy")
        normal = symbols("nx,ny")
        with pytest.raises(AssertionError):
            LineThroughPoint(point, direction=direction, normal=normal)

    def test_direction_specified(self):
        point = symbols("x,y")
        direction = symbols("dx,dy")
        line = LineThroughPoint(point, direction=direction)
        assert (line.dot(PointToVec3(point))).is_zero
        assert IdealPointOnLine(line) == Matrix(direction + (0,))

    def test_normal_vector_specified(self):
        point = symbols("x,y")
        normal = symbols("nx,ny")
        line = LineThroughPoint(point, normal=normal)
        assert (line.dot(PointToVec3(point))).is_zero
        assert line[0:2] == list(normal)


class TestPerpendicularBisector:
    def test_symbolic(self):
        x1, y1, x2, y2 = symbols("x1,y1,x2,y2")
        point1 = (x1, y1)
        point2 = (x2, y2)
        bisector = PerpendicularBisector(point1, point2)
        center = Matrix([(x1 + x2) / 2, (y1 + y2) / 2, 1])
        assert bisector.dot(center).expand().is_zero
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
