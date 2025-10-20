from sympy import Function, I, Matrix, nan, sqrt, symbols

from lib.matrix import (
    ConicMatrix,
    IsDefinite,
    IsNonZeroMultiple,
    IsPositiveMultiple,
    IsRealMatrix,
    MaxEigenvalue,
    MinEigenvalue,
    NonZeroCross,
    QuadraticForm,
)


class TestIsNonZeroMultiple:
    def test_vectors(self):
        assert IsNonZeroMultiple(Matrix([2]), Matrix([3])) is True
        assert IsNonZeroMultiple(Matrix([1, 2]), Matrix([-2, -4])) is True
        assert IsNonZeroMultiple(Matrix([1, 2]), Matrix([2, 1])) is False

    def test_zero_vectors(self):
        assert IsNonZeroMultiple(Matrix([0, 0]), Matrix([0, 0])) is True
        assert IsNonZeroMultiple(Matrix([1, 2]), Matrix([0, 0])) is False
        assert IsNonZeroMultiple(Matrix([0, 0]), Matrix([1, 2])) is False

    def test_matrices(self):
        matrix_1234 = Matrix([[1, 2], [3, 4]])
        assert IsNonZeroMultiple(Matrix([[1, 2]]), Matrix([[1], [2]])) is False
        assert IsNonZeroMultiple(matrix_1234, Matrix([[2, 4], [6, 8]])) is True
        assert IsNonZeroMultiple(matrix_1234, Matrix.zeros(2, 2)) is False
        assert IsNonZeroMultiple(matrix_1234, Matrix([[1, 1], [1, 1]])) is False

    def test_lists(self):
        assert IsNonZeroMultiple([1], [2, 3]) is False
        assert IsNonZeroMultiple([0, 0], [0, 0]) is True
        assert IsNonZeroMultiple([1, 2], [2, 4]) is True
        assert IsNonZeroMultiple([1, 2], [0, 0]) is False
        assert IsNonZeroMultiple([1, 2], [2, 3]) is False

    def test_matrix_vs_list(self):
        assert IsNonZeroMultiple(Matrix([1, 2]), [2, 4]) is True
        assert IsNonZeroMultiple([1, 2], Matrix([2, 4])) is True
        assert IsNonZeroMultiple(Matrix([1, 2]).T, [2, 4]) is False
        assert IsNonZeroMultiple(Matrix([1, 2]), [2, 3]) is False
        assert IsNonZeroMultiple([1, 2], Matrix([2, 3])) is False

    def test_symbolic(self):
        x = symbols("x")
        assert IsNonZeroMultiple([x, 1], [-x, -1]) is True

        p = symbols("p", positive=True)
        assert IsNonZeroMultiple([1, 2], [p, p * 2]) is True
        assert IsNonZeroMultiple([1, 2], [p, p]) is False

        nz = symbols("nz", nonzero=True)
        assert IsNonZeroMultiple([1, 2], [nz, nz * 2]) is True

    def test_undecidable(self):
        x = symbols("x")
        assert IsNonZeroMultiple([1, 2], [x, x * 2]) is None


class TestIsPositiveMultiple:
    def test_vectors(self):
        assert IsPositiveMultiple(Matrix([2]), Matrix([3])) is True
        assert IsPositiveMultiple(Matrix([1, 2]), Matrix([-2, -4])) is False
        assert IsPositiveMultiple(Matrix([1, 2]), Matrix([2, 1])) is False

    def test_zero_vectors(self):
        assert IsPositiveMultiple(Matrix([0, 0]), Matrix([0, 0])) is True
        assert IsPositiveMultiple(Matrix([1, 2]), Matrix([0, 0])) is False
        assert IsPositiveMultiple(Matrix([0, 0]), Matrix([1, 2])) is False

    def test_matrices(self):
        matrix_1234 = Matrix([[1, 2], [3, 4]])
        assert IsPositiveMultiple(Matrix([[1, 2]]), Matrix([[1], [2]])) is False
        assert IsPositiveMultiple(matrix_1234, Matrix([[2, 4], [6, 8]])) is True
        assert IsPositiveMultiple(matrix_1234, Matrix.zeros(2, 2)) is False
        assert IsPositiveMultiple(matrix_1234, Matrix([[1, 1], [1, 1]])) is False

    def test_lists(self):
        assert IsPositiveMultiple([1], [2, 3]) is False
        assert IsPositiveMultiple([0, 0], [0, 0]) is True
        assert IsPositiveMultiple([1, 2], [2, 4]) is True
        assert IsPositiveMultiple([1, 2], [0, 0]) is False
        assert IsPositiveMultiple([1, 2], [2, 3]) is False

    def test_matrix_vs_list(self):
        assert IsPositiveMultiple(Matrix([1, 2]), [2, 4]) is True
        assert IsPositiveMultiple([1, 2], Matrix([2, 4])) is True
        assert IsPositiveMultiple(Matrix([1, 2]).T, [2, 4]) is False
        assert IsPositiveMultiple(Matrix([1, 2]), [2, 3]) is False
        assert IsPositiveMultiple([1, 2], Matrix([2, 3])) is False

    def test_symbolic(self):
        x = symbols("x", real=True)
        assert IsPositiveMultiple([x, 1], [-x, -1]) is False

        p = symbols("p", positive=True)
        assert IsPositiveMultiple([1, 2], [p, p * 2]) is True
        assert IsPositiveMultiple([1, 2], [p, p]) is False

    def test_undecidable(self):
        x = symbols("x", real=True)
        assert IsPositiveMultiple([1, 2], [x, x * 2]) is None


class TestEigenvalues:
    def test_eigenvalues(self):
        matrix = Matrix([[1, 2], [2, 3]])
        assert MaxEigenvalue(matrix) == 2 + sqrt(5)
        assert MinEigenvalue(matrix) == 2 - sqrt(5)


class TestConicMatrix:
    def test_conic_matrix(self):
        assert ConicMatrix(1, 2, 3, 4, 5, 6) == Matrix(
            [[1, 2, 4], [2, 3, 5], [4, 5, 6]],
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


class TestIsRealMatrix:
    def test_numeric(self):
        assert IsRealMatrix(Matrix([[1, 2], [3, 4]])) is True
        assert IsRealMatrix(Matrix([[1, 2], [3, I]])) is False

    def test_symbolic(self):
        assert IsRealMatrix(Matrix([symbols("x", real=True)])) is True
        assert IsRealMatrix(Matrix([symbols("x")])) is None
        assert IsRealMatrix(Matrix([symbols("x", nonzero=True) * I])) is False


class TestIsDefinite:
    def test_non_square(self):
        assert IsDefinite(Matrix([1, 2])) is False

    def test_non_symmetric(self):
        assert IsDefinite(Matrix([[1, 2], [3, 4]])) is False

    def test_1x1_numeric(self):
        assert IsDefinite(Matrix([1]))
        assert IsDefinite(Matrix([-1]))
        assert IsDefinite(Matrix([0])) is False

    def test_1x1_symbolic(self):
        assert IsDefinite(Matrix([symbols("x", positive=True)]))
        assert IsDefinite(Matrix([symbols("x")])) is None

    def test_2x2_numeric(self):
        assert IsDefinite(Matrix.eye(2))
        assert IsDefinite(-Matrix.eye(2))
        assert IsDefinite(Matrix([[2, 1], [1, 2]]))
        assert IsDefinite(Matrix([[1, 2], [2, 1]])) is False
        assert IsDefinite(Matrix([[2, 1], [1, 0]])) is False
        assert IsDefinite(Matrix([[-2, 1], [1, -2]]))

    def test_2x2_symbolic(self):
        x = symbols("x", positive=True)
        assert IsDefinite(Matrix([[x, 0], [0, x]]))
        assert IsDefinite(Matrix([[x, 1], [1, x]])) is None
        assert IsDefinite(Matrix([[x, 1], [1, 0]])) is False
        assert IsDefinite(Matrix([[1, x], [x, -1]])) is False

    def test_3x3_diagonal_symbolic(self):
        x, y, z = symbols("x y z", positive=True)
        assert IsDefinite(Matrix.diag([x, y, z]))
        assert IsDefinite(Matrix.diag([-x, -y, -z]))
        assert IsDefinite(Matrix.diag([x, y, -z])) is False
        assert IsDefinite(Matrix.diag([x, y, 0])) is False

        r = symbols("r", real=True)
        assert IsDefinite(Matrix.diag([r, 1, -1])) is False
        assert IsDefinite(Matrix.diag([r, -1, -1])) is None  # definite iff r<0

        n = symbols("n", negative=True)
        assert IsDefinite(Matrix.diag([n, -1, -1])) is True

    def test_3x3_symbolic(self):
        a, b, c, d, e, f = symbols("a b c d e f")
        assert IsDefinite(Matrix([[a, b, d], [b, c, e], [d, e, f]])) is None
        assert IsDefinite(Matrix([[a, b, d], [b, c, e], [d, e, 0]])) is False
        assert IsDefinite(Matrix([[0, b, d], [b, c, e], [d, e, f]])) is False
