from sympy import Matrix
from sympy.abc import x

from lib.line import AreParallel, HorizontalLine, IDEAL_LINE, ParallelLine, Y_AXIS


class TestAreParallel:
    def test_numeric(self):
        assert AreParallel(HorizontalLine(1), HorizontalLine(2))
        assert AreParallel(HorizontalLine(1), Y_AXIS) is False
        assert AreParallel(Y_AXIS, IDEAL_LINE)
        assert AreParallel(Y_AXIS, ParallelLine(Y_AXIS, (1, 1)))

    def test_symbolic(self):
        assert AreParallel(Matrix([x, x * 2, 0]), Matrix([x * 2, x * 4, 0]))
        assert AreParallel(Matrix([x, x + x, 0]), Matrix([x, x + 1, 0])) is None
