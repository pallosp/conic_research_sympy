from sympy import Expr, Function, Integer, Matrix, S
from sympy.core.logic import fuzzy_and, fuzzy_not, fuzzy_or

from lib.matrix import IsDefinite


class ConicNormFactor(Function):
    """When the conic matrix `C` is multiplied by this value (`±1`), it will
    have the following properties:

    - For non-degenerate conics, the conic equation will evaluate to a positive
      number at the focus point(s), i.e. `(fx fy 1)ᵀ C (fx fy 1) > 0`.
    - For point conics, the conic equation will evaluate to ≤0 at all finite
      `(x, y, 1)` points.
    - For line pair conics, there is no preferred representation: this value will
      always be 1.
    - May return an unevaluated `sympy.Function` for symbolic conics whose type
      or determinant sign cannot be determined.
    """

    is_odd = True

    @classmethod
    def eval(cls, conic: Matrix) -> int | None:
        """Internal implementation. Call `ConicNormFactor(conic)` directly."""
        det = conic.det()
        if det.is_positive:
            return 1
        if det.is_negative:
            return -1

        # degenerate conic
        if det.is_zero:
            is_point = IsPointConic(conic)

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

    def _eval_Abs(self) -> Integer:
        return S.One

    def _eval_power(self, exponent: Expr) -> Expr | None:
        if exponent.is_even:
            return S.One
        if exponent.is_odd:
            return self
        return None


def IsDegenerate(conic: Matrix) -> bool | None:
    """Tells whether the conic is degenerate.

    Degenerate conics consist of a single projective point or a pair of
    projective lines. The zero matrix is also considered degenerate.
    Returns None if undecidable.
    """
    return conic.det().is_zero


def IsNonDegenerate(conic: Matrix) -> bool | None:
    """Tells whether the conic is non-degenerate.

    Non-degenerate conics include real or imaginary ellipses, parabolas and
    hyperbolas. Returns None if undecidable.
    """
    return conic.det().is_nonzero


def IsFiniteConic(conic: Matrix) -> bool | None:
    """Tells whether all points on the conic are finite.

    Returns None if undecidable.
    """
    return conic[:2, :2].det().factor().is_positive


def IsImaginaryEllipse(conic: Matrix) -> bool | None:
    """Tells Whether the conic is an ellipse with a real center and imaginary
    radii.

    Returns None if undecidable.
    """
    return IsDefinite(conic)


def IsEllipse(conic: Matrix) -> bool | None:
    """Tells whether the conic is an ellipse with real radii.

    Returns None if undecidable.
    """
    return fuzzy_and(
        [
            IsNonDegenerate(conic),
            conic[:2, :2].det().is_positive,
            fuzzy_not(IsDefinite(conic)),
        ],
    )


def IsCircle(conic: Matrix) -> bool | None:
    """Tells whether the conic is a circle.

    Returns None if undecidable.
    """
    a, _, _, b, c, _, _, _, _ = conic
    return fuzzy_and(
        [
            (a - c).expand().is_zero,
            b.is_zero,
            (a * conic.det()).is_negative,
        ],
    )


def IsParabola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a parabola.

    Returns None if undecidable.
    """
    return fuzzy_and([IsNonDegenerate(conic), conic[:2, :2].det().is_zero])


def IsHyperbola(conic: Matrix) -> bool | None:
    """Tells whether the conic is a hyperbola.

    Returns None if undecidable.
    """
    return fuzzy_and([IsNonDegenerate(conic), conic[:2, :2].det().is_negative])


def IsCircular(conic: Matrix) -> bool | None:
    """Tells whether there is a single center point around which the conic is
    invariant under all rotations.

    Circles, imaginary circles, zero-radius circles have such circular symmetry.
    Double ideal lines are not considered circular. Returns None if undecidable.
    """
    return fuzzy_and(
        [
            conic[0].is_nonzero,
            conic[3].is_zero,
            (conic[0] - conic[4]).expand().is_zero,
        ],
    )


def IsLinePair(conic: Matrix) -> bool | None:
    """Tells whether the conic is the union of two projective lines.

    Returns None if undecidable.
    """
    a, _, _, b, c, _, d, e, f = conic
    return fuzzy_and(
        [
            IsDegenerate(conic),
            fuzzy_not(conic.is_zero_matrix),
            # If all sqrt subexpressions in SplitToLines' implementation
            # are nonnegative, the conic splits to real lines.
            (b * b - a * c).factor().is_nonnegative,
            (d * d - a * f).factor().is_nonnegative,
            (e * e - c * f).factor().is_nonnegative,
        ],
    )


def IsDoubleLine(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of two coincident projective lines.

    Returns None if undecidable.
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


def IsPointConic(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of a single projective point.

    Returns None if undecidable.

    A conic is a point conic iff it's degenerate and splits to two lines with
    complex coordinates.
    """
    a, _, _, b, c, _, d, e, f = conic
    return fuzzy_and(
        [
            IsDegenerate(conic),
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


def IsFinitePointConic(conic: Matrix) -> bool | None:
    """Tells whether the conic consists of a single finite (Euclidean) point.

    Returns None if undecidable.
    """
    return fuzzy_and([IsDegenerate(conic), IsFiniteConic(conic)])
