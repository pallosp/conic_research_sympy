from sympy import Function, Matrix, nan, sqrt, symbols

from lib.matrix import (
    ConicMatrix,
    IsScalarMultiple,
    MaxEigenvalue,
    MinEigenvalue,
    NonZeroCross,
    QuadraticForm,
)


class TestIsScalarMultiple:
    def test_vectors(self):
        assert IsScalarMultiple(Matrix([2]), Matrix([3]))
        assert IsScalarMultiple(Matrix([1, 2]), Matrix([-2, -4]))
        assert IsScalarMultiple(Matrix([1, 2]), Matrix([0, 0]))
        assert not IsScalarMultiple(Matrix([1, 2]), Matrix([2, 1]))

    def test_matrices(self):
        assert not IsScalarMultiple(Matrix([[1, 2]]), Matrix([[1], [2]]))
        assert IsScalarMultiple(Matrix([[1, 2], [3, 4]]), Matrix([[2, 4], [6, 8]]))
        assert IsScalarMultiple(Matrix([[1, 2], [3, 4]]), Matrix([[0, 0], [0, 0]]))
        assert not IsScalarMultiple(Matrix([[1, 2], [3, 4]]), Matrix([[1, 1], [1, 1]]))


class TestEigenvalues:
    def test_eigenvalues(self):
        matrix = Matrix([[1, 2], [2, 3]])
        assert MaxEigenvalue(matrix) == 2 + sqrt(5)
        assert MinEigenvalue(matrix) == 2 - sqrt(5)


class TestConicMatrix:
    def test_conic_matrix(self):
        assert ConicMatrix(1, 2, 3, 4, 5, 6) == Matrix(
            [[1, 2, 4], [2, 3, 5], [4, 5, 6]]
        )


class TestQuadraticForm:
    def test_quadratic_form(self):
        # https://www.wolframalpha.com/input?i={{1,2}}*{{1,2},{3,4}}*{{1},{2}}
        assert QuadraticForm(Matrix([[1, 2], [3, 4]]), Matrix([1, 2])) == 27


class TestNonZeroCross:
    def test_full_matrix(self):
        matrix = Matrix([[1, 2], [3, 4]])
        assert NonZeroCross(matrix) == (Matrix([1, 3]), Matrix([1, 2]).T)

    def test_sparse_matrix(self):
        matrix = Matrix([[0, 1], [2, 0]])
        assert NonZeroCross(matrix) == (Matrix([1, 0]), Matrix([0, 1]).T)

    def test_zero_matrix(self):
        assert NonZeroCross(Matrix.zeros(2, 2)) == nan

    def test_undecidable(self):
        x = symbols("x")
        # Not evaluated
        assert isinstance(NonZeroCross(Matrix([x])), Function)

    def test_symbols_in_elements(self):
        x = symbols("x", nonnegative=True)
        assert NonZeroCross(Matrix([[x - x, 0], [1, 0]]))[1] == Matrix([1, 0]).T
        assert NonZeroCross(Matrix([[0, 0], [x + 1, 0]]))[0] == Matrix([0, x + 1])
