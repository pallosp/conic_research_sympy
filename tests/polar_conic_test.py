import pytest
from sympy import Matrix, pi, simplify, symbols

from lib.central_conic import central_conic_vertices, conic_center
from lib.circle import UNIT_CIRCLE
from lib.conic_classes import is_hyperbola
from lib.ellipse import ellipse
from lib.hyperbola import UNIT_HYPERBOLA
from lib.incidence import are_collinear, conic_contains_point
from lib.line import line_through_point
from lib.matrix import conic_matrix, is_nonzero_multiple
from lib.point import point_to_xy
from lib.polar_conic import (
    POLAR_UNIT_CIRCLE,
    PolarOrigin,
    angle_at_point,
    conic_from_polar_matrix,
    curvature_sign_at_angle,
    ellipse_to_polar_matrix,
    hyperbola_to_polar_matrix,
    point_at_angle,
    tangent_at_angle,
)
from lib.transform import homography_from_samples


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


class TestCurvatureSignAtAngle:
    def test_circle(self):
        circle = Matrix.eye(3)
        assert curvature_sign_at_angle(circle, 0) == 1
        assert curvature_sign_at_angle(circle, pi / 2) == 1

    def test_reverse_circle(self):
        circle = Matrix.diag([1, -1, 1])
        assert curvature_sign_at_angle(circle, 0) == -1

    def test_hyperbola(self):
        circle = ((1, 0), (0, 1), (-1, 0), (0, -1))
        hyperbola = ((1, 0), (1, 1, 0), (-1, 0), (1, -1, 0))
        polar_hyperbola = homography_from_samples(circle, hyperbola)
        assert curvature_sign_at_angle(polar_hyperbola, 0) == -1
        assert curvature_sign_at_angle(polar_hyperbola, pi / 2) == 0
        assert curvature_sign_at_angle(polar_hyperbola, pi) == 1
        assert curvature_sign_at_angle(polar_hyperbola, -pi / 2) == 0


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


class TestEllipseToPolarMatrix:
    def test_numeric_ellipse(self):
        e = ellipse((1, 2), 3, 4, r1_direction=(5, 6))

        start_points = [PolarOrigin.HORIZONTAL, PolarOrigin.VERTICAL]
        polar_ellipses = {}

        for start in start_points:
            p = ellipse_to_polar_matrix(e, start=start)
            polar_ellipses[start] = p
            assert is_nonzero_multiple(conic_from_polar_matrix(p), e)
            assert are_collinear(
                [
                    point_at_angle(p, 0),
                    conic_center(e),
                    point_at_angle(p, pi),
                ]
            )
            assert point_at_angle(p, symbols("a"))[2] == 1

        start_point_h = point_at_angle(polar_ellipses[PolarOrigin.HORIZONTAL], 0)
        assert start_point_h[1] == 2

        start_point_v = point_at_angle(polar_ellipses[PolarOrigin.VERTICAL], 0)
        assert start_point_v[0] == 1

    def test_unsupported_polar_origin(self):
        e = ellipse((0, 0), 2, 1)
        with pytest.raises(ValueError, match="Unsupported PolarOrigin"):
            ellipse_to_polar_matrix(e, start=PolarOrigin.IDEAL_POINT)


class TestHyperbolaToPolarMatrix:
    def test_numeric_hyperbola(self):
        hyperbola = conic_matrix(1, 2, 3, 4, 5, 6)
        assert is_hyperbola(hyperbola)

        start_points = [PolarOrigin.VERTEX]
        polar_hyperbolas = {}

        for start in start_points:
            polar_hyperbola = hyperbola_to_polar_matrix(hyperbola, start=start)
            polar_hyperbolas[start] = polar_hyperbola
            polar_hyperbola = polar_hyperbola.applyfunc(simplify)

            rebuilt_hyperbola = conic_from_polar_matrix(polar_hyperbola)
            rebuilt_hyperbola = rebuilt_hyperbola.applyfunc(simplify)

            assert is_nonzero_multiple(rebuilt_hyperbola, hyperbola)

        start_point = point_at_angle(polar_hyperbolas[PolarOrigin.VERTEX], 0)
        vertex = central_conic_vertices(hyperbola)[0]
        assert point_to_xy(start_point) == vertex

    def test_unsupported_polar_origin(self):
        with pytest.raises(ValueError, match="Unsupported PolarOrigin"):
            hyperbola_to_polar_matrix(UNIT_HYPERBOLA, start=PolarOrigin.COVERTEX)
