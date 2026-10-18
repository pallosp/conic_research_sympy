from sympy import I, Matrix, sqrt, symbols

from lib.matrix import (
    conic_matrix,
    is_definite_matrix,
    is_full_rank,
    is_nonzero_multiple,
    is_positive_multiple,
    is_real_matrix,
    max_eigenvalue,
    min_eigenvalue,
    outer_product_root,
    pivot_column,
    pivot_column_and_row,
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


class TestPivotColumn:
    def test_numeric(self):
        assert pivot_column(Matrix([[1, 2], [3, 4]])) == Matrix([1, 3])
        assert pivot_column(Matrix([[0, 1], [0, 2]])) == Matrix([1, 2])
        assert pivot_column(Matrix([[0, 1], [2, 0]])) == Matrix([1, 0])
        assert pivot_column(Matrix.zeros(3, 2)) == Matrix.zeros(3, 1)

    def test_symbolic(self):
        pos = symbols("p", positive=True)
        assert pivot_column(Matrix.diag([0, pos])) == Matrix([0, pos])

        x = symbols("x")
        assert pivot_column(Matrix.diag([x, 1])) == Matrix([0, 1])

        x = symbols("x", nonnegative=True)
        assert pivot_column(Matrix([[x - x, 0], [x + 1, 0]])) == Matrix([0, x + 1])

    def test_complex(self):
        x = symbols("x")
        assert pivot_column(Matrix([[x, I], [1, 0]])) == Matrix([I, 0])

    def test_undecidable(self):
        x = symbols("x")
        column = pivot_column(Matrix([[x, 0, 0], [0, x, 0]]))
        assert column.is_MatrixExpr
        assert column.shape == (2, 1)
        assert column.T.shape == (1, 2)
        assert (Matrix([[1, 2]]) * column).shape == (1, 1)
        assert column.subs(x, 3) == Matrix([3, 0])


class TestPivotColumnAndRow:
    def test_numeric(self):
        column, row = pivot_column_and_row(Matrix([[0, 1, 2], [3, 4, 5]]))
        assert column == Matrix([1, 4])
        assert row == Matrix([0, 1, 2])
        assert pivot_column_and_row(Matrix([[0, 0], [1, 2]]))[1] == Matrix([1, 2])

    def test_zero_matrix(self):
        column, row = pivot_column_and_row(Matrix.zeros(2, 3))
        assert column == Matrix.zeros(2, 1)
        assert row == Matrix.zeros(3, 1)

    def test_symbolic(self):
        x = symbols("x", nonnegative=True)
        assert pivot_column_and_row(Matrix([[x - x, 0], [1, 0]]))[1] == Matrix([1, 0])

    def test_complex(self):
        x = symbols("x")
        _, row = pivot_column_and_row(Matrix([[x, 0], [0, 1 + I]]))
        assert row == Matrix([0, 1 + I])

    def test_same_pivot_as_pivot_column(self):
        x = symbols("x")
        matrix = Matrix([[x, 0, 1], [2, 3, 4]])
        column, row = pivot_column_and_row(matrix)
        assert column == pivot_column(matrix) == Matrix([1, 4])
        assert row == Matrix([x, 0, 1])

    def test_undecidable(self):
        x = symbols("x")
        matrix = Matrix([[x, 0, 0], [0, x, 0]])
        column, row = pivot_column_and_row(matrix)
        assert column == pivot_column(matrix)
        assert row.is_MatrixExpr
        assert row.shape == (3, 1)
        assert column.subs(x, 3) == Matrix([3, 0])
        assert row.subs(x, 3) == Matrix([3, 0, 0])


class TestOuterProductRoot:
    def test_real_vector(self):
        v = Matrix([2, 0, -3])
        root = outer_product_root(v * v.T)
        assert root * root.T == v * v.T

    def test_first_diagonal_element_zero(self):
        v = Matrix([0, 1, 2])
        assert outer_product_root(v * v.T) == v

    def test_imaginary_vector(self):
        v = Matrix([I, 2 * I, 0])
        root = outer_product_root(v * v.T)
        assert root * root.T == v * v.T

    def test_complex_vector(self):
        # The only non-zero diagonal element, 2i, is not real
        v = Matrix([1 + I, 0])
        root = outer_product_root(v * v.T)
        assert root.expand() == v

    def test_symbolic_vector(self):
        x, y = symbols("x y", positive=True)
        v = Matrix([x, y])
        assert outer_product_root(v * v.T) == v

    def test_zero_matrix(self):
        assert outer_product_root(Matrix.zeros(3, 3)) == Matrix.zeros(3, 1)

    def test_undecidable(self):
        x = symbols("x")
        v = Matrix([x, 0])
        root = outer_product_root(v * v.T)
        assert root.is_MatrixExpr
        assert root.shape == (2, 1)
        assert str(root).startswith("OuterProductRoot(")
        assert root.subs(x, 3) == Matrix([3, 0])


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


class TestIsFullRank:
    def test_numeric_square_matrix(self):
        assert is_full_rank(Matrix([[1, 2], [3, 4]])) is True
        assert is_full_rank(Matrix([[1, 2], [2, 4]])) is False

    def test_numeric_rect_matrix(self):
        assert is_full_rank(Matrix([[0, 1, 2], [1, 2, 3]])) is True
        assert is_full_rank(Matrix([[0, 1, 2], [1, 2, 3]]).T) is True
        assert is_full_rank(Matrix([[1, 2, 3], [2, 4, 6]])) is False
        assert is_full_rank(Matrix([[1, 2, 3], [2, 4, 6]]).T) is False

    def test_symbolic_square_matrix(self):
        x = symbols("x")
        assert is_full_rank(Matrix([[x, x + 1], [x + 1, x + 2]])) is True
        assert is_full_rank(Matrix([[x, 2 * x], [3, 6]])) is False

    def test_symbolic_rect_matrix(self):
        x = symbols("x")

        full_rank_matrix = Matrix([[x, x + 1, x + 2], [x + 1, x + 2, x + 3]])
        assert is_full_rank(full_rank_matrix) is True
        assert is_full_rank(full_rank_matrix.T) is True

        lin_dep_matrix = Matrix([[x, x * 2, x * 3], [2, 4, 6]])
        assert is_full_rank(lin_dep_matrix) is False
        assert is_full_rank(lin_dep_matrix.T) is False

    def test_custom_simplifier(self):
        x, y = symbols("x y", real=True)
        m = Matrix([[1, sqrt(x)], [sqrt(y), sqrt(x * y)]])

        # The determinant can only be simplified to 0 when x ≥ 0 and y ≥ 0
        assert is_full_rank(m) is None
        assert is_full_rank(m, simplifier=lambda e: e.powsimp(force=True)) is False
