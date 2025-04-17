from sympy import simplify, symbols

from lib.central_conic import AxisLengths, ConicCenter
from lib.circle import UnitCircle
from lib.ellipse import Ellipse


class TestEllipseFromParams:
    def test_unit_circle_ellipse(self):
        assert Ellipse(0, 0, 1, 1) == UnitCircle()

    def test_center(self):
        x, y, r1, r2, angle = symbols("x,y,r1,r2,angle")
        ellipse = Ellipse(x, y, r1, r2, angle)
        center_x, center_y = ConicCenter(ellipse)
        assert x == simplify(center_x)
        assert y == simplify(center_y)

    def test_axis_lengths(self):
        x, y, angle = symbols("x,y,angle")
        r_min, r_diff = symbols("r_min,r_diff", nonnegative=True)
        ellipse = Ellipse(x, y, r_min, r_min + r_diff, angle)
        axes = [simplify(len) for len in AxisLengths(ellipse)]
        assert [r_min, r_min + r_diff] == axes or [r_min + r_diff, r_min] == axes
