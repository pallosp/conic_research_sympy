from sympy import Matrix, nan, sqrt, symbols

from lib.matrix import (
    ConicMatrix,
    IsScalarMultiple,
    MaxEigenvalue,
    MinEigenvalue,
    NonZeroCol,
    NonZeroRow,
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


class TestNonZeroRowOrColumn:
    def test_non_zero_row_numeric(self):
        assert NonZeroRow(Matrix([[1, 2], [3, 4]])) == Matrix([1, 2]).T
        assert NonZeroRow(Matrix([[0, 0], [0, 1]])) == Matrix([0, 1]).T
        assert NonZeroRow(Matrix([[0, 0], [0, 0]])) == nan

    def test_non_zero_column_numeric(self):
        assert NonZeroCol(Matrix([[1, 2], [3, 4]])) == Matrix([1, 3])
        assert NonZeroCol(Matrix([[0, 1], [0, 2]])) == Matrix([1, 2])
        assert NonZeroCol(Matrix([[0, 0], [0, 0]])) == nan

    def test_non_zero_row_symbolic(self):
        x = symbols("x")
        assert NonZeroRow(Matrix([[0, 0], [x, 0]])) == Matrix([x, 0]).T
        assert NonZeroRow(Matrix([[x - x, 0], [1, 0]])) == Matrix([1, 0]).T
