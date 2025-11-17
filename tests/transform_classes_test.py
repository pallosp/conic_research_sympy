from sympy import Matrix, expand, simplify, symbols

from lib.transform import reflect_to_line, rotate, translate
from lib.transform_classes import is_affine_transform, is_homography, is_similarity


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

    def test_simplifier(self):
        t = rotate(symbols("theta"))
        assert is_homography(t) is True
        assert is_homography(t, simplifier=expand) is None
        assert is_homography(t, simplifier=simplify) is True


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

    def test_simplifier(self):
        t = rotate(symbols("theta"))
        assert is_affine_transform(t) is True
        assert is_affine_transform(t, simplifier=expand) is None
        assert is_affine_transform(t, simplifier=simplify) is True


class TestIsSimilarity:
    def test_general_matrix(self):
        t = Matrix(3, 3, symbols("a b c d e f g h i"))
        assert is_similarity(t) is None

    def test_non_3x3_matrix(self):
        t = Matrix(2, 2, symbols("a b c d"))
        assert is_similarity(t) is False

    def test_affine_not_similarity(self):
        # A shear transformation is affine but not similarity
        t = Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
        assert is_affine_transform(t) is True
        assert is_similarity(t) is False

    def test_translate(self):
        t = translate(1, 2)
        assert is_similarity(t) is True

    def test_rotate(self):
        t = rotate(symbols("theta"))
        assert is_similarity(t) is True

    def test_scale(self):
        # Uniform scaling
        t = Matrix([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
        assert is_similarity(t) is True

        # Non-uniform scaling (not similarity)
        t = Matrix([[2, 0, 0], [0, 3, 0], [0, 0, 1]])
        assert is_similarity(t) is False

    def test_rotate_and_translate(self):
        t = rotate(0.5) * translate(1, 2)
        assert is_similarity(t) is True

    def test_scale_and_translate(self):
        t = Matrix([[2, 0, 1], [0, 2, 2], [0, 0, 1]])
        assert is_similarity(t) is True

    def test_symbolic_similarity(self):
        s, x, y = symbols("s x y")
        t = Matrix([[s, 0, x], [0, s, y], [0, 0, 1]])
        assert is_similarity(t) is None

    def test_symbolic_rotation(self):
        theta, x, y = symbols("theta x y")
        t = rotate(theta, x, y)
        assert is_similarity(t) is True

    def test_reflection(self):
        axis = Matrix(symbols("a b c", positive=True))
        t = reflect_to_line(axis)
        assert is_similarity(t) is True

    def test_simplifier(self):
        t = rotate(symbols("theta"))
        assert is_similarity(t) is True
        assert is_similarity(t, simplifier=expand) is None
        assert is_similarity(t, simplifier=simplify) is True
