from sympy import sqrt, symbols

from lib.conic import ConicFromPoly
from lib.hyperbola import HyperbolaFromFociAndPoint
from lib.matrix import IsNonZeroMultiple


class TestHyperbolaFromFociAndPoint:
    def test_numeric(self):
        # https://www.wolframalpha.com/input?i=x*y=4+focus
        f1 = (sqrt(8), sqrt(8))
        f2 = (-sqrt(8), -sqrt(8))
        p = (4, 1)
        hyperbola = HyperbolaFromFociAndPoint(f1, f2, p)
        x, y = symbols("x y")
        assert IsNonZeroMultiple(hyperbola, ConicFromPoly(x * y - 4))
