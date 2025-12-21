from typing import override

from sympy import Expr, Function, I, Integer, Matrix, S, sign, sqrt

from lib.conic_classes import is_point_conic


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

    def _eval_power(self, exponent: Expr) -> Expr | None:
        if exponent.is_even:
            return S.One
        if exponent.is_odd:
            return self
        return None


def focal_axis_direction(conic: Matrix) -> Matrix:
    """Returns the ideal point representing the direction of a conic's focal axis.

    Properties:
    - The focal axis is treated as an undirected line; its angle to the
      horizontal lies in the (-π/2, π/2] interval. For the full direction of a
      parabola, use [parabola_direction](#parabola.parabola_direction) instead.
    - Returns `[0, 0, 0]ᵀ` for circles and
      [circular conics](#conic_classification.is_circular).
    - Point conics constructed by
      [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(ellipse)
      preserve the axis direction of the original real or imaginary ellipse.
    - Line pair conics constructed by
      [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(hyperbola)
      have no such property.
    - The focal axis of `line_pair(l1, l2)`, and `angle_bisector(l1, l2)` point
      to the same direction.

    *Formula*:
    [research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)
    """
    # Use b as in ax² + 2bxy + cy²
    a, b, c = conic[0], conic[3], conic[4]
    norm_sign = ConicNormFactor(conic)
    x, y = sqrt(norm_sign * (a - c + 2 * I * b)).simplify().as_real_imag()
    return Matrix([x, y, 0])
