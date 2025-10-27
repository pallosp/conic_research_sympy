from collections.abc import Callable
from typing import override

from sympy import Expr, Function, Integer, Matrix, S, factor, sign
from sympy.core.logic import fuzzy_and, fuzzy_not, fuzzy_or

from lib.matrix import is_definite_matrix


class ConicNormFactor(Function):
    """Computes a normalization factor (±1) for a conic matrix `C`.

    When `C` is multiplied by this factor, the resulting conic has the
    following properties:

    - *Non-degenerate conics*: the conic equation evaluates to a positive
      value at the focus point(s), i.e. `[fx fy 1]ᵀ C [fx fy 1] > 0`.
    - *Point conics*: The conic equation evaluates to ≤0 for all finite
      points `[x, y, 1]ᵀ`.
    - *Line-pair conics*: no preferred normalization exists; the factor is
      always 1.
    - *Symbolic conics*: may return an unevaluated `sympy.Function` if the
      conic type or determinant sign cannot be determined.
    - `conic.det() * ConicNormFactor(conic) == Abs(conic.det())` holds for
      all conic types.
    """

    # Implies is_real = True, is_integer = True, and is_nonzero = True.
    is_odd = True

    @classmethod
    def eval(cls, conic: Matrix) -> int | None:  # noqa: PLR0911
        """Internal implementation. Call `ConicNormFactor(conic)` directly."""
        det = conic.det().factor()
        if det.is_positive:
            return 1
        if det.is_negative:
            return -1
        if det.is_real and det.is_nonzero:
            return sign(det)

        # degenerate conic
        if det.is_zero:
            is_point = is_point_conic(conic)

            # line pair
            if is_point is False:
                return 1

            # point conic
            if is_point is True:
                diag = conic.diagonal()
                if any(e.is_positive for e in diag):
                    return -1
                if any(e.is_negative for e in diag):
                    return 1

        return None

    @override
    def _eval_Abs(self) -> Integer:
        return S.One

    @override
    def _eval_power(self, exponent: Expr) -> Expr | None:
        if exponent.is_even:
            return S.One
        if exponent.is_odd:
            return self
        return None


def is_degenerate(
    conic: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = lambda expr: expr,
) -> bool | None:
    """Tells whether the conic is degenerate.

    Degenerate conics consist of a single projective point or a pair of
    projective lines. The zero matrix is also considered degenerate.

    Takes an optional `simplifier` callback that simplifies the conic matrix
    determinant before it gets compared to zero. Returns `None` if the result
    is undecidable.
    """
    return simplifier(conic.det()).is_zero


def is_nondegenerate(
    conic: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = lambda expr: expr,
) -> bool | None:
    """Tells whether the conic is non-degenerate.

    Non-degenerate conics include real or imaginary ellipses, parabolas and
    hyperbolas.

    Takes an optional `simplifier` callback that simplifies the conic matrix
    determinant before it gets compared to zero. Returns `None` if the result
    is undecidable.
    """
    return simplifier(conic.det()).is_nonzero


def is_central_conic(
    conic: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = lambda expr: expr,
) -> bool | None:
    """Tells whether a conic has a finite center of symmetry.

    Takes an optional `simplifier` callback that simplifies the central
    conicness polynomial before it gets compared to zero. Returns `None` if
    the result is undecidable.
    """
    return simplifier(conic[:2, :2].det()).is_nonzero


def is_finite_conic(
    conic: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = factor,
) -> bool | None:
    """Tells whether all points on the conic are finite.

    Takes an optional `simplifier` callback that simplifies the finiteness
    polynomial before it gets compared to zero. Returns `None` if the result
    is undecidable.
    """
    return simplifier(conic[:2, :2].det()).is_positive


def is_imaginary_ellipse(conic: Matrix) -> bool | None:
    """Tells whether the conic is an imaginary ellipse.

    Imaginary ellipses have real center and focus points, but imaginary radii
    and eccentricity. In addition, all solutions of their conic equations are
    points with complex coordinates.

    Returns `None` if undecidable.
    """
    return is_definite_matrix(conic)


def is_ellipse(conic: Matrix) -> bool | None:
    """Tells whether the conic is an ellipse.

    Returns `False` for
    [imaginary ellipses](#conic_classification.is_imaginary_ellipse),
    or `None` if the conic's type is undecidable.
    """
    return fuzzy_and(
        [
            is_nondegenerate(conic),
            conic[:2, :2].det().is_positive,
            fuzzy_not(is_definite_matrix(conic)),
        ],
    )


def is_circle(conic: Matrix) -> bool | None:
    """Tells whether the conic is a circle.

    Returns `None` if undecidable.
    """
    a, _, _, b, c, _, _, _, _ = conic
    return fuzzy_and(
        [
            (a - c).expand().is_zero,
            b.is_zero,
            (a * conic.det()).is_negative,
        ],
    )


def is_parabola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a parabola.

    Returns `None` if undecidable.
    """
    return fuzzy_and([is_nondegenerate(conic), conic[:2, :2].det().is_zero])


def is_hyperbola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a hyperbola.

    Returns `None` if undecidable.
    """
    return fuzzy_and([is_nondegenerate(conic), conic[:2, :2].det().is_negative])


def is_circular(conic: Matrix) -> bool | None:
    """Tells whether there is a single center point around which the conic is
    invariant under all rotations.

    Circles, imaginary circles, zero-radius circles have such circular symmetry.
    Double ideal lines are not considered circular. Returns `None` if undecidable.
    """
    return fuzzy_and(
        [
            conic[0].is_nonzero,
            conic[3].is_zero,
            (conic[0] - conic[4]).expand().is_zero,
        ],
    )


def is_line_pair(conic: Matrix) -> bool | None:
    """Tells whether the conic is the union of two projective lines.

    Returns `None` if undecidable.
    """
    a, _, _, b, c, _, d, e, f = conic
    return fuzzy_and(
        [
            is_degenerate(conic),
            fuzzy_not(conic.is_zero_matrix),
            # If all sqrt subexpressions in SplitToLines' implementation
            # are nonnegative, the conic splits to real lines.
            (b * b - a * c).factor().is_nonnegative,
            (d * d - a * f).factor().is_nonnegative,
            (e * e - c * f).factor().is_nonnegative,
        ],
    )


def is_double_line(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of two coincident projective lines.

    Returns `None` if undecidable.
    """
    # The conic represents a double line iff its matrix has rank 1.
    #
    # A matrix has rank ≤ 1 iff the cross product of every pair of columns is
    # zero, or equivalently its adjugate is the zero matrix.
    return fuzzy_and(
        [
            conic.adjugate().is_zero_matrix,
            fuzzy_not(conic.is_zero_matrix),
        ],
    )


def is_point_conic(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of a single projective point.

    Returns `None` if undecidable.

    A conic is a point conic iff it's degenerate and splits to two lines with
    complex coordinates.
    """
    a, _, _, b, c, _, d, e, f = conic
    return fuzzy_and(
        [
            is_degenerate(conic),
            # If any of the sqrt subexpressions in SplitToLines' implementation
            # are negative, the conic splits to complex lines.
            fuzzy_or(
                [
                    (b * b - a * c).factor().is_negative,
                    (d * d - a * f).factor().is_negative,
                    (e * e - c * f).factor().is_negative,
                ],
            ),
        ],
    )


def is_finite_point_conic(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of a single finite (Euclidean) point.

    Returns `None` if undecidable.
    """
    return fuzzy_and([is_degenerate(conic), is_finite_conic(conic)])
