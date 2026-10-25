from collections.abc import Callable, Sequence
from typing import Self, cast, overload

from sympy import (
    Expr,
    Function,
    ImmutableMatrix,
    Integer,
    Matrix,
    MatrixBase,
    MatrixExpr,
    S,
    ShapeError,
    expand,
    sqrt,
)
from sympy.core.logic import fuzzy_and, fuzzy_not
from sympy.matrices.expressions.matexpr import MatrixElement
from sympy.printing.pretty.pretty import PrettyPrinter
from sympy.printing.pretty.stringpict import prettyForm

from lib.type_aliases import Conic, Scalar


def is_nonzero_multiple(
    m1: MatrixBase | MatrixExpr | Sequence[Scalar],
    m2: MatrixBase | MatrixExpr | Sequence[Scalar],
) -> bool | None:
    """Tells whether two matrices are non-zero scalar multiples of each other.

    Treats lists and tuples as column vectors, and compares unevaluated matrix
    expressions element-wise. Returns `None` if undecidable.
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


def is_positive_multiple(
    m1: MatrixBase | Sequence[Scalar],
    m2: MatrixBase | Sequence[Scalar],
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


def max_eigenvalue(symmetric_matrix_2x2: MatrixBase) -> Expr:
    """Returns the higher eigenvalue of a 2x2 symmetric matrix."""
    m = symmetric_matrix_2x2
    if m.shape != (2, 2) or m != m.T:
        raise ValueError("The input must be a 2x2 symmetric matrix.")
    a, b, _, c = m.flat()
    return (a + c) / 2 + sqrt((a - c) ** 2 + 4 * b**2) / 2


def min_eigenvalue(symmetric_matrix_2x2: MatrixBase) -> Expr:
    """Returns the lower eigenvalue of a 2x2 symmetric matrix."""
    m = symmetric_matrix_2x2
    if m.shape != (2, 2) or m != m.T:
        raise ValueError("The input must be a 2x2 symmetric matrix.")
    a, b, _, c = m.flat()
    return (a + c) / 2 - sqrt((a - c) ** 2 + 4 * b**2) / 2


def conic_matrix(  # noqa: PLR0913, PLR0917 (too-many-arguments, too-many-positional-arguments)
    a: Scalar,
    b: Scalar,
    c: Scalar,
    d: Scalar,
    e: Scalar,
    f: Scalar,
) -> Conic:
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


def quadratic_form(sym_matrix: MatrixBase, vector: MatrixBase) -> Expr:
    """Computes the quadratic form for a n⨯n symmetric matrix and an n-element
    column vector.

    Use case: when `sym_matrix` and `vector` represent a conic and a projective
    point, respectively, the quadratic form is zero iff the point is on the
    conic.

    *Formula*: `vᵀ·M·v`
    """
    return vector.dot(vector.T * sym_matrix)


def _pivot_index(matrix: MatrixBase) -> int | None:
    """Finds the first provably non-zero element of a matrix.

    Returns 0 for zero matrices, or `None` if undecidable.
    """
    all_zero = True
    for i, el in enumerate(matrix.flat()):
        is_zero = el.is_zero
        if is_zero is False:
            return i
        if is_zero is None:
            all_zero = False
    return 0 if all_zero else None


class MatrixFunction(MatrixExpr):
    """Base class of matrix valued functions of a matrix, which stay
    unevaluated until the result can be determined.

    Subclasses implement
     - `eval(cls, matrix)` as a class method, returning the result or `None` if
       it depends on undecidable assumptions, like `sympy.Function.eval`;
     - the `shape` property of the result.

    The argument may be an unevaluated `MatrixExpr` itself. Unevaluated
    instances are `MatrixExpr`s too. Substitutions such as `subs` rebuild them,
    which retries the evaluation.
    """

    eval: Callable[[MatrixBase | MatrixExpr], MatrixBase | None]

    def __new__(
        cls,
        matrix: MatrixBase | MatrixExpr,
        *,
        evaluate: bool = True,
    ) -> "MatrixBase | Self":
        """Evaluates the function if possible, unless `evaluate` is false."""
        if evaluate:
            result = cls.eval(matrix)
            if result is not None:
                return result
        if isinstance(matrix, MatrixBase):
            matrix = ImmutableMatrix(matrix)
        return super().__new__(cls, matrix)

    @property
    def matrix(self) -> MatrixBase | MatrixExpr:
        """The argument of the function."""
        return cast("MatrixBase | MatrixExpr", self.args[0])

    def _entry(self, i: Expr, j: Expr, **_kwargs: object) -> Expr:
        return MatrixElement(self, i, j)

    def _pretty(self, printer: PrettyPrinter) -> prettyForm:
        # Without this, the pretty printer would fall back to str().
        return printer._print_Function(  # noqa: SLF001 (private-member-access)
            Function(type(self).__name__)(*self.args)
        )


class SkewMatrix(MatrixFunction):
    """Unevaluated form of [skew_matrix](#matrix.skew_matrix)."""

    @classmethod
    def eval(cls, vector3: MatrixBase | MatrixExpr) -> MatrixBase | None:
        """Internal implementation. Call `skew_matrix(vector3)` instead."""
        if vector3.shape != (3, 1):
            msg = f"Expected a 3D column vector, got shape {vector3.shape}"
            raise ShapeError(msg)
        if not isinstance(vector3, MatrixBase):
            return None
        return vector3.hat()

    @property
    def shape(self) -> tuple[Integer, Integer]:
        """The shape of the result: `(3, 3)`."""
        return (Integer(3), Integer(3))

    def _entry(self, i: Expr, j: Expr, **_kwargs: object) -> Expr:
        if not (i.is_Integer and j.is_Integer):
            return MatrixElement(self, i, j)
        if i == j:
            return S.Zero
        # Same layout as Matrix.hat()
        element = self.matrix[int(3 - i - j)]
        return element if (i - j) % 3 == 1 else -element


@overload
def skew_matrix(vector3: MatrixBase) -> MatrixBase: ...
@overload
def skew_matrix(vector3: MatrixExpr) -> SkewMatrix: ...
def skew_matrix(vector3: MatrixBase | MatrixExpr) -> MatrixBase | SkewMatrix:
    """Creates a skew-symmetric matrix from a 3D column vector.

    Returns an explicit matrix for an explicit vector, and an unevaluated
    `SkewMatrix` for a matrix expression.
    """
    return SkewMatrix(vector3)


class PivotColumn(MatrixFunction):
    """Unevaluated form of [pivot_column](#matrix.pivot_column)."""

    @classmethod
    def eval(cls, matrix: MatrixBase | MatrixExpr) -> MatrixBase | None:
        """Internal implementation. Call `pivot_column(matrix)` instead."""
        matrix = matrix.as_explicit()
        index = _pivot_index(matrix)
        return None if index is None else matrix.col(index % matrix.cols)

    @property
    def shape(self) -> tuple[Integer, Integer]:
        """The shape of the result: `(matrix.rows, 1)`."""
        return (Integer(self.matrix.rows), S.One)


class TransposedPivotRow(MatrixFunction):
    """Unevaluated form of the transposed row returned by
    [pivot_column_and_row](#matrix.pivot_column_and_row).
    """

    @classmethod
    def eval(cls, matrix: MatrixBase | MatrixExpr) -> MatrixBase | None:
        """Internal implementation. Call `pivot_column_and_row(matrix)` instead."""
        matrix = matrix.as_explicit()
        index = _pivot_index(matrix)
        return None if index is None else matrix.row(index // matrix.cols).T

    @property
    def shape(self) -> tuple[Integer, Integer]:
        """The shape of the result: `(matrix.cols, 1)`."""
        return (Integer(self.matrix.cols), S.One)


def pivot_column(matrix: MatrixBase | MatrixExpr) -> MatrixBase | PivotColumn:
    """Returns the column of a matrix for its first provably nonzero element.

    The result is a zero vector in case of a zero matrix.

    Returns an unevaluated column vector shaped `MatrixExpr` if none of the
    elements can be proven non-zero.

    [pivot_column_and_row](#matrix.pivot_column_and_row) also returns the row
    of the same pivot.
    """
    return PivotColumn(matrix)


def pivot_column_and_row(
    matrix: MatrixBase | MatrixExpr,
) -> tuple[MatrixBase | PivotColumn, MatrixBase | TransposedPivotRow]:
    """Returns the column and the transposed row of a matrix for its first
    provably nonzero element, i.e. two column vectors.

    The column is the same as [pivot_column](#matrix.pivot_column)`(matrix)`.

    Returns unevaluated column vector shaped `MatrixExpr`s if none of the
    elements can be proven non-zero.
    """
    explicit = matrix.as_explicit()
    index = _pivot_index(explicit)
    if index is not None:
        return (
            explicit.col(index % explicit.cols),
            explicit.row(index // explicit.cols).T,
        )
    return (
        PivotColumn(matrix, evaluate=False),
        TransposedPivotRow(matrix, evaluate=False),
    )


class OuterProductRoot(MatrixFunction):
    """Unevaluated form of [outer_product_root](#matrix.outer_product_root)."""

    @classmethod
    def eval(cls, matrix: MatrixBase | MatrixExpr) -> MatrixBase | None:
        """Internal implementation. Call `outer_product_root(matrix)` instead."""
        matrix = matrix.as_explicit()
        diagonal = matrix.diagonal()
        index = _pivot_index(diagonal)
        if index is None:
            return None
        if diagonal[index].is_zero:
            # All diagonal elements are zero, therefore v = 0.
            return matrix.col(0)
        return matrix.col(index) / sqrt(diagonal[index])

    @property
    def shape(self) -> tuple[Integer, Integer]:
        """The shape of the result: `(matrix.rows, 1)`."""
        return (Integer(self.matrix.rows), S.One)


def outer_product_root(
    matrix: MatrixBase | MatrixExpr,
) -> MatrixBase | OuterProductRoot:
    """Finds a vector `v` such that `v·vᵀ` equals the given symmetric matrix.

    Such a vector exists iff the matrix has rank ≤ 1. It is unique up to sign,
    and may be complex. The result is unspecified for other matrices.

    Returns a zero vector in case of a zero matrix, or an unevaluated column
    vector shaped `MatrixExpr` if none of the diagonal elements can be proven
    to be non-zero.

    *Formula*: `M.col(i) / √Mᵢᵢ` for any `i` with `Mᵢᵢ ≠ 0`, since
    `(v·vᵀ).col(i) = v·vᵢ` and `(v·vᵀ)ᵢᵢ = vᵢ²`.
    """
    return OuterProductRoot(matrix)


def is_real_matrix(matrix: MatrixBase) -> bool | None:
    """Checks if all elements of a matrix are real.

    Returns a bool or `None` if undecidable.
    """
    return fuzzy_and(el.is_real for el in matrix.flat())


def is_definite_matrix(matrix: MatrixBase) -> bool | None:
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
    diag = matrix.diagonal().flat()
    has_ge_0_diag_el = any(e.is_nonnegative for e in diag)
    has_le_0_diag_el = any(e.is_nonpositive for e in diag)
    if has_ge_0_diag_el and has_le_0_diag_el:
        return False

    # A matrix is positive definite iff all leading principal minors > 0
    size = matrix.rows
    return fuzzy_and(matrix[:i, :i].det().is_positive for i in range(1, size + 1))


def is_full_rank(
    matrix: MatrixBase,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether a matrix has full rank.

    Takes an optional `simplifier` callback that simplifies the determinant
    before it gets compared to zero. Returns `None` if the result is
    undecidable.
    """
    rows, cols = matrix.shape
    if all(el.is_Rational for el in matrix.flat()):
        return matrix.rank() == min(rows, cols)

    if cols < rows:
        matrix = matrix.T * matrix
    elif cols > rows:
        matrix = matrix * matrix.T

    # The Laplace method is usually faster for symbolic matrices than the
    # default Bareiss algorithm.
    return simplifier(matrix.det("laplace")).is_nonzero
