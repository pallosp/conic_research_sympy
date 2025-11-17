from sympy import Matrix, symbols

from lib.transform import rotate, translate
from lib.transform_classes import is_affine_transform, is_homography


class TestIsHomography:
    def test_general_matrix(self):
        t = Matrix(3, 3, symbols("a b c d e f g h i"))
        assert is_homography(t) is None

    def test_singular_matrix(self):
        t = Matrix([[1, 2, 3], [4, 5, 6], [5, 7, 9]])  # Row 3 = Row 1 + Row 2
        assert t.det() == 0
        assert is_homography(t) is False

    def test_non_singular_matrix(self):
        t = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        assert t.det() != 0
        assert is_homography(t) is True

    def test_non_3x3_matrix(self):
        t = Matrix(2, 2, symbols("a b c d"))
        assert is_homography(t) is False


class TestIsAffineTransform:

    def test_general_matrix(self):
        t = Matrix(3, 3, symbols("a b c d e f g h i"))
        assert is_affine_transform(t) is None

    def test_potentially_singular(self):
        t = Matrix(3, 3, [*symbols("a b c d e f"), 0, 0, 1])
        assert is_affine_transform(t) is None

    def test_translate(self):
        t = translate(*symbols("x y"))
        assert is_affine_transform(t) is True

    def test_arbitrary_bottom_right_element(self):
        t = translate(*symbols("x y")) * -2
        assert is_affine_transform(t) is True

    def test_rotate(self):
        t = rotate(symbols("a"), symbols("x"), symbols("y"))
        assert is_affine_transform(t) is True
