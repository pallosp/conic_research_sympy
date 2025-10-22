from collections.abc import Sequence

from sympy import Expr, Function, Matrix, nan, sqrt
from sympy.core.logic import fuzzy_and, fuzzy_not
from sympy.core.numbers import NaN


def IsNonZeroMultiple(
    m1: Matrix | Sequence[Expr],
    m2: Matrix | Sequence[Expr],
) -> bool | None:
    """Tells whether two matrices are non-zero scalar multiples of each other.

    Treats lists and tuples as column vectors. Returns `None` if undecidable.
    """
    if not isinstance(m1, Matrix):
        m1 = Matrix(m1)
    if not isinstance(m2, Matrix):
        m2 = Matrix(m2)

    if m1.shape != m2.shape:
        return False

    m1_is_zero = m1.is_zero_matrix
    m2_is_zero = m2.is_zero_matrix
    if m1_is_zero and m2_is_zero:
        return True

    # Linearize the matrices to be able to use dot product.
    v1 = m1.reshape(len(m1), 1)
    v2 = m2.reshape(len(m2), 1)

    return fuzzy_and(
        [
            (v1.dot(v1) * v2.dot(v2) - v1.dot(v2) ** 2).simplify().is_zero,
            fuzzy_not(m1_is_zero),
            fuzzy_not(m2_is_zero),
        ],
    )


def IsPositiveMultiple(
    m1: Matrix | Sequence[Expr],
    m2: Matrix | Sequence[Expr],
) -> bool | None:
    """Tells whether two matrices are positive scalar multiples of each other.

    Treats lists and tuples as column vectors. Returns `None` if undecidable.
    """
    if not isinstance(m1, Matrix):
        m1 = Matrix(m1)
    if not isinstance(m2, Matrix):
        m2 = Matrix(m2)

    if m1.shape != m2.shape:
        return False

    m1_is_zero = m1.is_zero_matrix
    m2_is_zero = m2.is_zero_matrix
    if m1_is_zero and m2_is_zero:
        return True

    # Linearize the matrices to be able to use dot product.
    v1 = m1.reshape(len(m1), 1)
    v2 = m2.reshape(len(m2), 1)

    return fuzzy_and(
        [
            v1.dot(v2).factor().is_positive,
            (v1.dot(v1) * v2.dot(v2) - v1.dot(v2) ** 2).simplify().is_zero,
        ],
    )


def MaxEigenvalue(symmetric_matrix_2x2: Matrix) -> Expr:
    """Returns the higher eigenvalue of a 2x2 symmetric matrix."""
    m = symmetric_matrix_2x2
    if m.shape != (2, 2) or m != m.T:
        raise ValueError("The input must be a 2x2 symmetric matrix.")
    a, b, _, c = m
    return (a + c) / 2 + sqrt((a - c) ** 2 + 4 * b**2) / 2


def MinEigenvalue(symmetric_matrix_2x2: Matrix) -> Expr:
    """Returns the lower eigenvalue of a 2x2 symmetric matrix."""
    m = symmetric_matrix_2x2
    if m.shape != (2, 2) or m != m.T:
        raise ValueError("The input must be a 2x2 symmetric matrix.")
    a, b, _, c = m
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
    """Computes the quadratic form for a n⨯n symmetric matrix and an n-element
    column vector.

    Use case: when `sym_matrix` and `vector` represent a conic and a projective
    point, respectively, the quadratic form is zero iff the point is on the
    conic.

    *Formula*: `vᵀ·M·v`
    """
    return vector.dot(vector.T * sym_matrix)


def SkewMatrix(vector3: Matrix) -> Matrix:
    """Creates a skew-symmetric matrix from a 3D vector."""
    x, y, z = vector3
    return Matrix(
        [
            [0, -z, y],
            [z, 0, -x],
            [-y, x, 0],
        ],
    )


class NonZeroCross(Function):
    """Finds a column and a row in a matrix whose intersection is a non-zero
    element.

    Returns an unevaluated `sympy.Function` if none of the elements can be
    proven to be non-zero, or `nan` in case of a zero matrix.
    """

    @classmethod
    def eval(cls, matrix: Matrix) -> tuple[Matrix, Matrix] | NaN | None:
        """Internal implementation. Call `NonZeroCross(matrix)` directly."""
        all_zero = True
        for i in range(len(matrix)):
            el = matrix[i]
            if el.is_nonzero:
                return (matrix.col(i % matrix.cols), matrix.row(i // matrix.cols))
            if all_zero and not el.is_zero:
                all_zero = False
        if all_zero:
            return nan
        return None


def IsRealMatrix(matrix: Matrix) -> bool | None:
    """Checks if all elements of a matrix are real.

    Returns a bool or `None` if undecidable.
    """
    return fuzzy_and(el.is_real for el in matrix)


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
