from sympy import Matrix, symbols

from lib.central_conic import ConicCenter, SemiMajorAxis, SemiMinorAxis
from lib.circle import Circle


class TestConicCenter:
    def test_circle(self):
        x, y, r = symbols("x,y,r")
        circle = Circle(x, y, r)
        center_x, center_y = ConicCenter(circle)
        assert x == center_x
        assert y == center_y


class TestAxes:
    def test_circle_radius(self):
        x, y = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = Circle(x, y, r)
        assert r == SemiMajorAxis(circle)
        assert r == SemiMinorAxis(circle)
