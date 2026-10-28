import itertools
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

    return _is_nonzero_multiple(m1.flat(), m2.flat())


def _is_nonzero_multiple(a: list[Expr], b: list[Expr]) -> bool | None:
    """Tells whether two vectors of the same length are non-zero scalar multiples
    of each other.
    """
    a_is_zero = fuzzy_and(e.is_zero for e in a)
    b_is_zero = fuzzy_and(e.is_zero for e in b)
    if a_is_zero and b_is_zero:
        return True
    if a_is_zero or b_is_zero:
        return False

    # The vectors are proportional iff all 2x2 minors of [a | b] vanish.
    # Unlike the Lagrange identity |a|²|b|² - (a·b)² = Σ minor², this also
    # works over complex numbers.
    def is_zero_minor(i: int, j: int) -> bool | None:
        return _is_zero(a[i] * b[j] - a[j] * b[i])

    # If a[k] or b[k] ≠ 0, the minors containing index k imply proportionality.
    pivot = next(
        (k for k in range(len(a)) if a[k].is_zero is False or b[k].is_zero is False),
        None,
    )
    if pivot is not None:
        # fuzzy_and stops at the first False.
        pivot_minors_vanish = fuzzy_and(
            is_zero_minor(pivot, j) for j in range(len(a)) if j != pivot
        )
        if pivot_minors_vanish is False:
            return False
        if pivot_minors_vanish:
            # a and b are proportional, so each is zero iff its pivot element is.
            return fuzzy_and(
                [fuzzy_not(_is_zero(a[pivot])), fuzzy_not(_is_zero(b[pivot]))],
            )

    # The other minors can't prove proportionality, but may disprove it.
    pairs = itertools.combinations(range(len(a)), 2)
    if any(is_zero_minor(i, j) is False for i, j in pairs if pivot not in (i, j)):
        return False
    return None


def _is_zero(expr: Expr) -> bool | None:
    """Checks whether an expression is zero, trying the cheap `expand` first."""
    is_zero = expr.expand().is_zero
    return is_zero if is_zero is not None else expr.simplify().is_zero


def _is_positive(expr: Expr) -> bool | None:
    """Checks whether an expression is positive, factoring it if necessary."""
    is_positive = expr.is_positive
    return is_positive if is_positive is not None else expr.factor().is_positive


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

    if m1.is_zero_matrix and m2.is_zero_matrix:
        return True

    a = m1.flat()
    b = m2.flat()

    # If b = λ·a, then λ = b[k] / a[k] for any a[k] ≠ 0, and 1/λ = a[k] / b[k]
    # for any b[k] ≠ 0. Both have the same sign, and checking it is cheap, so do
    # it before the proportionality check. The ratios are equal, but SymPy may
    # only be able to decide the sign for some of them.
    ratio_is_positive = None
    for ak, bk in zip(a, b, strict=True):
        if ak.is_zero is False:
            ratio_is_positive = _is_positive(bk / ak)
        elif bk.is_zero is False:
            ratio_is_positive = _is_positive(ak / bk)
        if ratio_is_positive is not None:
            break
    if ratio_is_positive is False:
        return False

    return fuzzy_and([ratio_is_positive, _is_nonzero_multiple(a, b)])


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


def conic_coefficients(conic: Conic) -> tuple[Expr, Expr, Expr, Expr, Expr, Expr]:
    """Returns the elements `(a, b, c, d, e, f)` of a conic matrix.

    It's the inverse of [conic_matrix](#matrix.conic_matrix).
    """
    a, _, _, b, c, _, d, e, f = conic.flat()
    return a, b, c, d, e, f


def quadratic_form(
    sym_matrix: MatrixBase | MatrixExpr, vector: MatrixBase | MatrixExpr
) -> Expr:
    """Computes the quadratic form for a n⨯n symmetric matrix and an n-element
    column vector.

    Use case: when `sym_matrix` and `vector` represent a conic and a projective
    point, respectively, the quadratic form is zero iff the point is on the
    conic.

    Accepts unevaluated `MatrixExpr` arguments, whose elements appear as
    `MatrixElement`s in the result.

    *Formula*: `vᵀ·M·v`
    """
    # SymPy's stubs type MatrixExpr products as Expr (fix upstream).
    return (vector.T * sym_matrix * vector)[0, 0]  # pyright: ignore[reportIndexIssue]


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

    # mypy requires __new__ to return an instance of cls, but like
    # sympy.Function, this evaluates to a different type when possible.
    def __new__(  # type: ignore[misc]
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

    The matrix turns cross products into matrix products: `skew(v)·w = v⨯w`

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
