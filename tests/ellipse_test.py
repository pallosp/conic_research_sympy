from sympy import factor, pi, simplify, symbols

from lib.central_conic import AxisLengths, ConicCenter
from lib.circle import UNIT_CIRCLE
from lib.ellipse import Ellipse
from lib.matrix import IsScalarMultiple
from lib.point import ORIGIN


class TestEllipseFromParams:
    def test_unit_circle_ellipse(self):
        assert Ellipse(ORIGIN, 1, 1) == UNIT_CIRCLE

    def test_center(self):
        center = symbols("x,y")
        r1, r2 = symbols("r1,r2")
        r1_direction = symbols("r1_x,r1_y")
        ellipse = Ellipse(center, r1, r2, r1_direction=r1_direction)
        assert center == factor(ConicCenter(ellipse))

    def test_axis_direction_vector_length_invariance(self):
        center = symbols("x,y")
        r1, r2, r1_dir_x, r1_dir_y = symbols("r1,r2,r1_dir_x,r1_dir_y")
        conic1 = Ellipse(center, r1, r2, r1_direction=(r1_dir_x, r1_dir_y))
        conic2 = Ellipse(center, r1, r2, r1_direction=(r1_dir_x * -2, r1_dir_y * -2))
        assert IsScalarMultiple(conic1, conic2)

    def test_axis_lengths_circle(self):
        center = symbols("x,y")
        r = symbols("r", nonnegative=True)
        circle = Ellipse(center, r, r)
        assert AxisLengths(circle) == (r, r)

    def test_axis_lengths_numeric(self):
        ellipse = Ellipse((1, 2), 3, 4, r1_direction=(5, 6))
        assert AxisLengths(ellipse) == (3, 4)
        ellipse = Ellipse((1, 2), 3, 4, r1_angle=pi / 6)
        assert AxisLengths(ellipse) == (3, 4)

    def test_axis_lengths_general_case(self):
        center = symbols("x,y")
        r_min, r_diff = symbols("r_min,r_diff", positive=True)
        ellipse = Ellipse(center, r_min, r_min + r_diff, r1_direction=(73, -25))
        axes = tuple(factor(simplify(len)) for len in AxisLengths(ellipse))
        assert (r_min, r_min + r_diff) == axes or (r_min + r_diff, r_min) == axes
