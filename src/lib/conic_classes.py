from collections.abc import Callable

from sympy import Expr, Matrix, factor
from sympy.core.logic import fuzzy_and, fuzzy_not, fuzzy_or

from lib.matrix import is_definite_matrix


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
    [imaginary ellipses](#conic_classes.is_imaginary_ellipse),
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


def is_rectangular_hyperbola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a rectangular hyperbola.

    It's also called equilateral hyperbola or right hyperbola.

    *Definition and properties*:
    <https://mathworld.wolfram.com/RectangularHyperbola.html>

    Returns `None` if undecidable.
    """
    a, c, _ = conic.diagonal()
    return fuzzy_and([is_nondegenerate(conic), (a + c).is_zero])


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


#: The x² - y² = 1 hyperbola.
#: https://en.wikipedia.org/wiki/Unit_hyperbola
UNIT_HYPERBOLA: Matrix = Matrix.diag([1, -1, -1])
