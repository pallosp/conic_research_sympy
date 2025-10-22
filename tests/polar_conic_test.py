from sympy import Matrix, pi

from lib.circle import UNIT_CIRCLE
from lib.incidence import conic_contains_point
from lib.polar_conic import POLAR_UNIT_CIRCLE, conic_from_polar_matrix, point_at_angle


class TestConicFromPolarMatrix:
    def test_unit_circle(self):
        assert conic_from_polar_matrix(POLAR_UNIT_CIRCLE) == UNIT_CIRCLE

    def test_five_points(self):
        polar = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 0]])
        conic = conic_from_polar_matrix(polar)
        assert conic_contains_point(conic, point_at_angle(polar, 0))
        assert conic_contains_point(conic, point_at_angle(polar, pi / 4))
        assert conic_contains_point(conic, point_at_angle(polar, pi / 2))
        assert conic_contains_point(conic, point_at_angle(polar, pi))
        assert conic_contains_point(conic, point_at_angle(polar, -pi / 2))
