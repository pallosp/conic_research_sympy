from sympy import Matrix, Rational, symbols

from lib.polygon import signed_area


class TestSignedArea:
    def test_fewer_than_three_points(self):
        assert signed_area() == 0
        assert signed_area((1, 2)) == 0
        assert signed_area((1, 2), (3, 5)) == 0

    def test_counterclockwise_square(self):
        assert signed_area((0, 0), (2, 0), (2, 2), (0, 2)) == 4

    def test_clockwise_square(self):
        assert signed_area((0, 0), (0, 2), (2, 2), (2, 0)) == -4

    def test_homogeneous_coordinates(self):
        assert signed_area((0, 0, 1), (4, 0, 2), (3, 3, 3)) == 1

    def test_full_precision(self):
        assert signed_area((0, 0), (1, 0), (0, 1)) == Rational(1, 2)

    def test_symbolic_triangle(self):
        x1, y1, x2, y2, x3, y3 = symbols("x1 y1 x2 y2 x3 y3")
        det = Matrix([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]).det()
        assert signed_area((x1, y1), (x2, y2), (x3, y3)) * 2 == det

    def test_self_intersecting(self):
        # A bowtie whose loops have opposite orientations
        assert signed_area((0, 0), (2, 2), (2, 0), (0, 2)) == 0
