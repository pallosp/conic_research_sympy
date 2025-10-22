from sympy import Function, I, Matrix, nan, sqrt, symbols

from lib.matrix import (
    NonzeroCross,
    conic_matrix,
    is_definite_matrix,
    is_nonzero_multiple,
    is_positive_multiple,
    is_real_matrix,
    max_eigenvalue,
    min_eigenvalue,
    quadratic_form,
)


class TestIsNonZeroMultiple:
    def test_vectors(self):
        assert is_nonzero_multiple(Matrix([2]), Matrix([3])) is True
        assert is_nonzero_multiple(Matrix([1, 2]), Matrix([-2, -4])) is True
        assert is_nonzero_multiple(Matrix([1, 2]), Matrix([2, 1])) is False

    def test_zero_vectors(self):
        assert is_nonzero_multiple(Matrix([0, 0]), Matrix([0, 0])) is True
        assert is_nonzero_multiple(Matrix([1, 2]), Matrix([0, 0])) is False
        assert is_nonzero_multiple(Matrix([0, 0]), Matrix([1, 2])) is False

    def test_matrices(self):
        matrix_1234 = Matrix([[1, 2], [3, 4]])
        assert is_nonzero_multiple(Matrix([[1, 2]]), Matrix([[1], [2]])) is False
        assert is_nonzero_multiple(matrix_1234, Matrix([[2, 4], [6, 8]])) is True
        assert is_nonzero_multiple(matrix_1234, Matrix.zeros(2, 2)) is False
        assert is_nonzero_multiple(matrix_1234, Matrix([[1, 1], [1, 1]])) is False

    def test_lists(self):
        assert is_nonzero_multiple([1], [2, 3]) is False
        assert is_nonzero_multiple([0, 0], [0, 0]) is True
        assert is_nonzero_multiple([1, 2], [2, 4]) is True
        assert is_nonzero_multiple([1, 2], [0, 0]) is False
        assert is_nonzero_multiple([1, 2], [2, 3]) is False

    def test_matrix_vs_list(self):
        assert is_nonzero_multiple(Matrix([1, 2]), [2, 4]) is True
        assert is_nonzero_multiple([1, 2], Matrix([2, 4])) is True
        assert is_nonzero_multiple(Matrix([1, 2]).T, [2, 4]) is False
        assert is_nonzero_multiple(Matrix([1, 2]), [2, 3]) is False
        assert is_nonzero_multiple([1, 2], Matrix([2, 3])) is False

    def test_symbolic(self):
        x = symbols("x")
        assert is_nonzero_multiple([x, 1], [-x, -1]) is True

        p = symbols("p", positive=True)
        assert is_nonzero_multiple([1, 2], [p, p * 2]) is True
        assert is_nonzero_multiple([1, 2], [p, p]) is False

        nz = symbols("nz", nonzero=True)
        assert is_nonzero_multiple([1, 2], [nz, nz * 2]) is True

    def test_undecidable(self):
        x = symbols("x")
        assert is_nonzero_multiple([1, 2], [x, x * 2]) is None


class TestIsPositiveMultiple:
    def test_vectors(self):
        assert is_positive_multiple(Matrix([2]), Matrix([3])) is True
        assert is_positive_multiple(Matrix([1, 2]), Matrix([-2, -4])) is False
        assert is_positive_multiple(Matrix([1, 2]), Matrix([2, 1])) is False

    def test_zero_vectors(self):
        assert is_positive_multiple(Matrix([0, 0]), Matrix([0, 0])) is True
        assert is_positive_multiple(Matrix([1, 2]), Matrix([0, 0])) is False
        assert is_positive_multiple(Matrix([0, 0]), Matrix([1, 2])) is False

    def test_matrices(self):
        matrix_1234 = Matrix([[1, 2], [3, 4]])
        assert is_positive_multiple(Matrix([[1, 2]]), Matrix([[1], [2]])) is False
        assert is_positive_multiple(matrix_1234, Matrix([[2, 4], [6, 8]])) is True
        assert is_positive_multiple(matrix_1234, Matrix.zeros(2, 2)) is False
        assert is_positive_multiple(matrix_1234, Matrix([[1, 1], [1, 1]])) is False

    def test_lists(self):
        assert is_positive_multiple([1], [2, 3]) is False
        assert is_positive_multiple([0, 0], [0, 0]) is True
        assert is_positive_multiple([1, 2], [2, 4]) is True
        assert is_positive_multiple([1, 2], [0, 0]) is False
        assert is_positive_multiple([1, 2], [2, 3]) is False

    def test_matrix_vs_list(self):
        assert is_positive_multiple(Matrix([1, 2]), [2, 4]) is True
        assert is_positive_multiple([1, 2], Matrix([2, 4])) is True
        assert is_positive_multiple(Matrix([1, 2]).T, [2, 4]) is False
        assert is_positive_multiple(Matrix([1, 2]), [2, 3]) is False
        assert is_positive_multiple([1, 2], Matrix([2, 3])) is False

    def test_symbolic(self):
        x = symbols("x", real=True)
        assert is_positive_multiple([x, 1], [-x, -1]) is False

        p = symbols("p", positive=True)
        assert is_positive_multiple([1, 2], [p, p * 2]) is True
        assert is_positive_multiple([1, 2], [p, p]) is False

    def test_undecidable(self):
        x = symbols("x", real=True)
        assert is_positive_multiple([1, 2], [x, x * 2]) is None


class TestEigenvalues:
    def test_eigenvalues(self):
        matrix = Matrix([[1, 2], [2, 3]])
        assert max_eigenvalue(matrix) == 2 + sqrt(5)
        assert min_eigenvalue(matrix) == 2 - sqrt(5)


class TestConicMatrix:
    def test_conic_matrix(self):
        assert conic_matrix(1, 2, 3, 4, 5, 6) == Matrix(
            [[1, 2, 4], [2, 3, 5], [4, 5, 6]],
        )


class TestQuadraticForm:
    def test_quadratic_form(self):
        # https://www.wolframalpha.com/input?i={{1,2}}*{{1,2},{3,4}}*{{1},{2}}
        assert quadratic_form(Matrix([[1, 2], [3, 4]]), Matrix([1, 2])) == 27


class TestNonZeroCross:
    def test_full_matrix(self):
        matrix = Matrix([[1, 2], [3, 4]])
        assert NonzeroCross(matrix) == (Matrix([1, 3]), Matrix([1, 2]).T)

    def test_sparse_matrix(self):
        matrix = Matrix([[0, 1], [2, 0]])
        assert NonzeroCross(matrix) == (Matrix([1, 0]), Matrix([0, 1]).T)

    def test_zero_matrix(self):
        assert NonzeroCross(Matrix.zeros(2, 2)) == nan

    def test_undecidable(self):
        x = symbols("x")
        # Not evaluated
        assert isinstance(NonzeroCross(Matrix([x])), Function)

    def test_symbols_in_elements(self):
        x = symbols("x", nonnegative=True)
        assert NonzeroCross(Matrix([[x - x, 0], [1, 0]]))[1] == Matrix([1, 0]).T
        assert NonzeroCross(Matrix([[0, 0], [x + 1, 0]]))[0] == Matrix([0, x + 1])


class TestIsRealMatrix:
    def test_numeric(self):
        assert is_real_matrix(Matrix([[1, 2], [3, 4]])) is True
        assert is_real_matrix(Matrix([[1, 2], [3, I]])) is False

    def test_symbolic(self):
        assert is_real_matrix(Matrix([symbols("x", real=True)])) is True
        assert is_real_matrix(Matrix([symbols("x")])) is None
        assert is_real_matrix(Matrix([symbols("x", nonzero=True) * I])) is False


class TestIsDefinite:
    def test_non_square(self):
        assert is_definite_matrix(Matrix([1, 2])) is False

    def test_non_symmetric(self):
        assert is_definite_matrix(Matrix([[1, 2], [3, 4]])) is False

    def test_1x1_numeric(self):
        assert is_definite_matrix(Matrix([1]))
        assert is_definite_matrix(Matrix([-1]))
        assert is_definite_matrix(Matrix([0])) is False

    def test_1x1_symbolic(self):
        assert is_definite_matrix(Matrix([symbols("x", positive=True)]))
        assert is_definite_matrix(Matrix([symbols("x")])) is None

    def test_2x2_numeric(self):
        assert is_definite_matrix(Matrix.eye(2))
        assert is_definite_matrix(-Matrix.eye(2))
        assert is_definite_matrix(Matrix([[2, 1], [1, 2]]))
        assert is_definite_matrix(Matrix([[1, 2], [2, 1]])) is False
        assert is_definite_matrix(Matrix([[2, 1], [1, 0]])) is False
        assert is_definite_matrix(Matrix([[-2, 1], [1, -2]]))

    def test_2x2_symbolic(self):
        x = symbols("x", positive=True)
        assert is_definite_matrix(Matrix([[x, 0], [0, x]]))
        assert is_definite_matrix(Matrix([[x, 1], [1, x]])) is None
        assert is_definite_matrix(Matrix([[x, 1], [1, 0]])) is False
        assert is_definite_matrix(Matrix([[1, x], [x, -1]])) is False

    def test_3x3_diagonal_symbolic(self):
        x, y, z = symbols("x y z", positive=True)
        assert is_definite_matrix(Matrix.diag([x, y, z]))
        assert is_definite_matrix(Matrix.diag([-x, -y, -z]))
        assert is_definite_matrix(Matrix.diag([x, y, -z])) is False
        assert is_definite_matrix(Matrix.diag([x, y, 0])) is False

        r = symbols("r", real=True)
        assert is_definite_matrix(Matrix.diag([r, 1, -1])) is False
        assert is_definite_matrix(Matrix.diag([r, -1, -1])) is None  # definite iff r<0

        n = symbols("n", negative=True)
        assert is_definite_matrix(Matrix.diag([n, -1, -1])) is True

    def test_3x3_symbolic(self):
        a, b, c, d, e, f = symbols("a b c d e f")
        assert is_definite_matrix(Matrix([[a, b, d], [b, c, e], [d, e, f]])) is None
        assert is_definite_matrix(Matrix([[a, b, d], [b, c, e], [d, e, 0]])) is False
        assert is_definite_matrix(Matrix([[0, b, d], [b, c, e], [d, e, f]])) is False
