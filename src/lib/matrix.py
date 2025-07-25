from typing import Sequence, Tuple
from sympy import Expr, Function, Matrix, nan, sqrt
from sympy.core.logic import fuzzy_and
from sympy.core.numbers import NaN


def IsNonZeroMultiple(m1: Matrix | Sequence[Expr], m2: Matrix | Sequence[Expr]) -> bool:
    """Tells whether two matrices are non-zero scalar multiples of each other.

    Treats lists and tuples as column vectors.
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


def MaxEigenvalue(symmetric_matrix2x2: Matrix) -> Expr:
    """Returns the higher eigenvalue of a 2x2 symmetric matrix."""
    assert symmetric_matrix2x2.shape == (2, 2)
    a, b, b2, c = symmetric_matrix2x2
    assert b == b2
    return (a + c) / 2 + sqrt((a - c) ** 2 + 4 * b**2) / 2


def MinEigenvalue(symmetric_matrix2x2: Matrix):
    """Returns the lower eigenvalue of a 2x2 symmetric matrix."""
    assert symmetric_matrix2x2.shape == (2, 2)
    a, b, b2, c = symmetric_matrix2x2
    assert b == b2
    return (a + c) / 2 - sqrt((a - c) ** 2 + 4 * b**2) / 2


def ConicMatrix(a: Expr, b: Expr, c: Expr, d: Expr, e: Expr, f: Expr) -> Matrix:
    """Builds a 3x3 symmetric conic matrix from its elements.

    The conic equation looks like this:
    ```
            [a b d] [x]
    [x y 1] [b c e] [y] = 0
            [d e f] [1]
    ```

    or in expanded form `ax² + 2bxy + cy² + 2dx + 2ey + f = 0`
    """
    return Matrix([[a, b, d], [b, c, e], [d, e, f]])


def QuadraticForm(sym_matrix: Matrix, vector: Matrix) -> Expr:
    """Computes the quadratic form for a nxn symmetric matrix and an n-element
    column vector.

    Use case: when `sym_matrix` and `vector` represent a conic and a projective
    point, respectively, the quadratic form is zero iff the point is on the
    conic.

    Formula: `vᵀ·M·v`
    """
    return vector.dot(vector.T * sym_matrix)


def SkewMatrix(vector3: Matrix) -> Matrix:
    """Creates a skew-symmetric matrix from a 3d vector."""
    x, y, z = vector3
    return Matrix(
        [
            [0, -z, y],
            [z, 0, -x],
            [-y, x, 0],
        ]
    )


class NonZeroCross(Function):
    """Finds a column and a row in a matrix whose intersection is a non-zero
    element.

    Returns an unevaluated `sympy.Function` if none of the elements can be
    proven to be non-zero, or `nan` in case of a zero matrix.
    """

    @classmethod
    def eval(cls, matrix: Matrix) -> Tuple[Matrix, Matrix] | NaN | None:
        all_zero = True
        for i in range(len(matrix)):
            el = matrix[i]
            if el.is_nonzero:
                return (matrix.col(i % matrix.cols), matrix.row(i // matrix.cols))
            if all_zero and not el.is_zero:
                all_zero = False
        if all_zero:
            return nan


def IsDefinite(matrix: Matrix) -> bool | None:
    """Checks if a real matrix is either positive or negative definite.

    Returns a bool or `None` if the definitess can't be decided.
    """
    # The simple alternative algorithm,
    #
    #   fuzzy_or([matrix.is_positive_definite, matrix.is_negative_definite])
    #
    # would be able to prove definiteness or non-definiteness of symbolic
    # matrices in fewer cases.

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
