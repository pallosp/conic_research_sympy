from sympy import Matrix, symbols
from sympy.abc import x, y

from lib.line import (
    AreParallel,
    ArePerpendicular,
    HorizontalLine,
    IDEAL_LINE,
    LineBetween,
    ParallelLine,
    PerpendicularBisector,
    PerpendicularLine,
    VerticalLine,
    X_AXIS,
    Y_AXIS,
)


class TestPerpendicularLine:
    def test_numeric(self):
        assert Y_AXIS == PerpendicularLine(X_AXIS, (0, 1))


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
