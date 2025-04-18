from sympy import symbols

from lib.central_conic import ConicCenter, SemiMajorAxis, SemiMinorAxis
from lib.circle import Circle, UNIT_CIRCLE
from lib.transform import ScaleXY, TransformConic


class TestConicCenter:
    def test_circle(self):
        x, y, r = symbols("x,y,r")
        circle = Circle((x, y), r)
        assert (x, y) == ConicCenter(circle)


class TestAxes:
    def test_circle_radius(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = Circle(center, r)
        assert r == SemiMajorAxis(circle)
        assert r == SemiMinorAxis(circle)

    def test_ellipse_axes(self):
        ellipse = TransformConic(UNIT_CIRCLE, ScaleXY(2, 3))
        assert SemiMajorAxis(ellipse) == 3
        assert SemiMajorAxis(ellipse * -1) == 3
        assert SemiMinorAxis(ellipse) == 2
        assert SemiMinorAxis(ellipse * -1) == 2
