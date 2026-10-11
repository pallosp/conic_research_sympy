from lib.central_conic import conic_center
from lib.incidence import conic_contains_point
from lib.point import centroid, point_to_vec3
from lib.triangle import steiner_ellipse, steiner_inellipse


class TestSteinerEllipse:
    def test_circumellipse(self):
        p1 = (4, 1)
        p2 = (7, 3)
        p3 = (5, 5)
        ellipse = steiner_ellipse(p1, p2, p3)
        assert conic_contains_point(ellipse, point_to_vec3(p1))
        assert conic_contains_point(ellipse, point_to_vec3(p2))
        assert conic_contains_point(ellipse, point_to_vec3(p3))
        assert conic_center(ellipse) == centroid(p1, p2, p3)

    def test_inellipse(self):
        p1 = (4, 1)
        p2 = (7, 3)
        p3 = (5, 5)
        ellipse = steiner_inellipse(p1, p2, p3)
        assert conic_contains_point(ellipse, point_to_vec3(centroid(p2, p3)))
        assert conic_contains_point(ellipse, point_to_vec3(centroid(p3, p1)))
        assert conic_contains_point(ellipse, point_to_vec3(centroid(p1, p2)))
        assert conic_center(ellipse) == centroid(p1, p2, p3)
