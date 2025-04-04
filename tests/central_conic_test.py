from sympy import Matrix
from sympy import symbols

from lib.central_conic import ConicCenter
from lib.circle import Circle


class TestConicCenter:
    def test_circle(self):
        x, y, r = symbols("x,y,r")
        circle = Circle(x, y, r)
        center_x, center_y = ConicCenter(circle)
        assert x == center_x
        assert y == center_y
