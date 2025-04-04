from sympy import Matrix, sqrt

from lib.matrix import IsScalarMultiple, MaxEigenvalue, MinEigenvalue


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
