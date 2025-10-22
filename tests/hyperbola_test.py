from sympy import sqrt, symbols

from lib.conic import conic_from_poly
from lib.hyperbola import hyperbola_from_foci_and_point
from lib.matrix import is_nonzero_multiple


class TestHyperbolaFromFociAndPoint:
    def test_numeric(self):
        # https://www.wolframalpha.com/input?i=x*y=4+focus
        f1 = (sqrt(8), sqrt(8))
        f2 = (-sqrt(8), -sqrt(8))
        p = (4, 1)
        hyperbola = hyperbola_from_foci_and_point(f1, f2, p)
        x, y = symbols("x y")
        assert is_nonzero_multiple(hyperbola, conic_from_poly(x * y - 4))
