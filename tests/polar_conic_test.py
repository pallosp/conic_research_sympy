from sympy import Matrix, pi

from lib.circle import UNIT_CIRCLE
from lib.incidence import conic_contains_point
from lib.line import line_through_point
from lib.matrix import is_nonzero_multiple
from lib.polar_conic import (
    POLAR_UNIT_CIRCLE,
    angle_at_point,
    conic_from_polar_matrix,
    point_at_angle,
    tangent_at_angle,
)


class TestAngleAtPoint:
    polar_conic = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 0]])
    for angle in [0, pi / 4, pi / 2, pi, -pi / 2]:
        point = point_at_angle(polar_conic, angle)
        assert angle == angle_at_point(polar_conic, point)


class TestTangentAtAngle:
    def test_unit_circle(self):
        tangent = tangent_at_angle(POLAR_UNIT_CIRCLE, pi / 4)
        point = point_at_angle(POLAR_UNIT_CIRCLE, pi / 4)
        expected = line_through_point(point, direction=(-1, 1))
        assert is_nonzero_multiple(tangent, expected)


class TestConicFromPolarMatrix:
    def test_unit_circle(self):
        assert conic_from_polar_matrix(POLAR_UNIT_CIRCLE) == UNIT_CIRCLE

    def test_five_points(self):
        polar_conic = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 0]])
        conic = conic_from_polar_matrix(polar_conic)
        assert conic_contains_point(conic, point_at_angle(polar_conic, 0))
        assert conic_contains_point(conic, point_at_angle(polar_conic, pi / 4))
        assert conic_contains_point(conic, point_at_angle(polar_conic, pi / 2))
        assert conic_contains_point(conic, point_at_angle(polar_conic, pi))
        assert conic_contains_point(conic, point_at_angle(polar_conic, -pi / 2))
