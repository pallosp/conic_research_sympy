from sympy import Matrix

from lib.intersection import ConicXLine
from lib.matrix import ConicMatrix, QuadraticForm


class TestConicXLine:
    def test_conic_x_line(self):
        conic = ConicMatrix(1, 2, 3, 4, 5, 6)
        line = Matrix([1, 2, 3])
        p1, p2 = ConicXLine(conic, line)
        assert QuadraticForm(conic, p1).equals(0)
        assert QuadraticForm(conic, p2).equals(0)
