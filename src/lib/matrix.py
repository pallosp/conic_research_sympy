from sympy import Function, Matrix, nan, sqrt
from sympy.core.logic import fuzzy_and


def IsNonZeroMultiple(m1: Matrix | list, m2: Matrix | list) -> bool:
    """Tells whether two matrices are non-zero scalar multiples of each other.

    Treats lists as column vectors.
    """
    if not isinstance(m1, Matrix):
        m1 = Matrix(m1)
    if not isinstance(m2, Matrix):
        m2 = Matrix(m2)
    if m1.shape != m2.shape:
        return False
    if m1.is_zero_matrix != m2.is_zero_matrix:
        return False
    v1 = m1.reshape(len(m1), 1)
    v2 = m2.reshape(len(m2), 1)
    return (v1.dot(v1) * v2.dot(v2) - v1.dot(v2) ** 2).equals(0)


def MaxEigenvalue(symmetric_matrix2x2: Matrix):
    assert symmetric_matrix2x2.shape == (2, 2)
    a, b, b2, c = symmetric_matrix2x2
    assert b == b2
    return (a + c) / 2 + sqrt((a - c) ** 2 + 4 * b**2) / 2


def MinEigenvalue(symmetric_matrix2x2: Matrix):
    assert symmetric_matrix2x2.shape == (2, 2)
    a, b, b2, c = symmetric_matrix2x2
    assert b == b2
    return (a + c) / 2 - sqrt((a - c) ** 2 + 4 * b**2) / 2


def ConicMatrix(a, b, c, d, e, f):
    """3x3 symmetric matrix from the conic equation

            [a b d] [x]
    [x y 1] [b c e] [y] = 0
            [d e f] [1]

    Expanded form: ax² + 2bxy + cy² + 2dx + 2ey + f = 0
    """
    return Matrix([[a, b, d], [b, c, e], [d, e, f]])


def QuadraticForm(sym_matrix: Matrix, vector: Matrix):
    """Quadratic form for a nxn symmetric matrix and a n-element column vector.

    Formula: vᵀ·M·v
    """
    return vector.dot(vector.T * sym_matrix)


def SkewMatrix(vector3: Matrix):
    """Skew-symmetric matrix for a 3d vector."""
    x, y, z = vector3
    return Matrix(
        [
            [0, -z, y],
            [z, 0, -x],
            [-y, x, 0],
        ]
    )


class NonZeroCross(Function):
    @classmethod
    def eval(cls, matrix: Matrix):
        all_zero = True
        for i in range(len(matrix)):
            el = matrix[i]
            if el.is_nonzero:
                return (matrix.col(i % matrix.cols), matrix.row(i // matrix.cols))
            if all_zero and not el.is_zero:
                all_zero = False
        if all_zero:
            return nan


def IsDefinite(matrix: Matrix) -> bool:
    """Checks if a real matrix is either positive or negative definite."""
    # Non-square or non-symmetric matrices are not definite
    if not matrix.is_symmetric():
        return False

    if matrix[0].is_negative:
        matrix = -matrix

    # Take a quick look at the diagonal to potentially prove non-definiteness
    has_ge_0_diag_el = any(e.is_nonnegative for e in matrix.diagonal())
    has_le_0_diag_el = any(e.is_nonpositive for e in matrix.diagonal())
    if has_ge_0_diag_el and has_le_0_diag_el:
        return False

    # A matrix is positive definite iff all leading principal minors > 0
    size = matrix.rows
    return fuzzy_and(matrix[:i, :i].det().is_positive for i in range(1, size + 1))
