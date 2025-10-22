from sympy import Matrix, pi

from lib.circle import UNIT_CIRCLE
from lib.incidence import ConicContainsPoint
from lib.polar_conic import POLAR_UNIT_CIRCLE, ConicFromPolarMatrix, PointAtAngle


class TestConicFromPolarMatrix:
    def test_unit_circle(self):
        assert ConicFromPolarMatrix(POLAR_UNIT_CIRCLE) == UNIT_CIRCLE

    def test_five_points(self):
        polar = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 0]])
        conic = ConicFromPolarMatrix(polar)
        assert ConicContainsPoint(conic, PointAtAngle(polar, 0))
        assert ConicContainsPoint(conic, PointAtAngle(polar, pi / 4))
        assert ConicContainsPoint(conic, PointAtAngle(polar, pi / 2))
        assert ConicContainsPoint(conic, PointAtAngle(polar, pi))
        assert ConicContainsPoint(conic, PointAtAngle(polar, -pi / 2))
