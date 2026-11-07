from sympy import (
    Expr,
    I,
    Integer,
    Piecewise,
    S,
    sign,
    sqrt,
)
from sympy.core.function import DefinedFunction

from lib.matrix import conic_coefficients
from lib.point import ideal_point
from lib.type_aliases import Conic, PointVector


class ConicNormFactor(DefinedFunction):
    """Computes a normalization factor (±1) for a conic matrix `C`.

    When `C` is multiplied by this factor, the resulting conic has the
    following properties:

    - *Non-degenerate conics*: the conic equation evaluates to a positive
      value at the focus point(s), i.e. `[fx fy 1]ᵀ C [fx fy 1] > 0`. For
      parabolas, this makes the upper left 2×2 block negative semidefinite
      ([proof](../src/research/conic_properties/normalized_parabola.py)).
    - *Point conics*: the matrix is negative semidefinite, i.e. the conic
      equation evaluates to ≤0 at all points.
    - *Parallel line pairs and double lines*: like parabolas, these conics
      have a zero discriminant. Therefore the upper left 2×2 block of their
      matrix is singular, hence semidefinite. The factor makes it negative
      semidefinite, as it is for normalized parabolas.

      This keeps the normalization continuous: as `ε → 0`, the normalized
      matrices of the parabolas `y² = εx + k` converge to those of parallel
      lines if `k > 0`, and to that of a double line if `k = 0`.
    - *Crossing line pairs, conics containing the ideal line*: the factor
      is 1, since no preferred normalization exists.
    - *Symbolic conics*: may remain unevaluated if undecidable.

    `C.det() * ConicNormFactor(C) == Abs(C.det())` holds for all conic types.
    """

    # Implies is_real = True, is_integer = True, and is_nonzero = True.
    is_odd = True

    @classmethod
    def eval(cls, conic: Conic) -> Expr | None:
        """Internal implementation. Call `ConicNormFactor(conic)` directly."""
        det = conic.det().factor()
        if det.is_nonzero:  # implies det.is_real
            return sign(det)  # ±1 if det's sign is known

        if det.is_zero:
            return _degenerate_conic_norm_factor(conic)

        # It's undecidable whether the conic is degenerate. The result is
        # valid if sign(det) and the degenerate conic's factor agree.
        if det.is_nonnegative and _degenerate_conic_norm_factor(conic) == 1:
            return S.One

        if det.is_nonpositive and _degenerate_conic_norm_factor(conic) == -1:
            return S.NegativeOne

        return None

    def _eval_Abs(self) -> Integer:  # noqa: N802 (hook name looked up by Abs)
        return S.One

    def _eval_power(self, exponent: Expr) -> Expr | None:
        if exponent.is_even:
            return S.One
        if exponent.is_odd:
            return self
        return None


def _degenerate_conic_norm_factor(conic: Conic) -> Expr | None:
    """Computes `ConicNormFactor` for degenerate conics.

    If the discriminant `b² - a·c` is ≤ 0, the upper left 2×2 block is
    semidefinite, and the factor makes it negative semidefinite. This covers
    point conics (making the whole matrix negative semidefinite), parallel
    line pairs and double lines. A semidefinite block has the sign of its
    trace.

    The factor is 1 for crossing line pairs (positive discriminant) and for
    conics containing the ideal line (zero block, zero trace).
    """
    a, b, c, *_ = conic_coefficients(conic)
    trace = a + c
    if trace.is_nonpositive:
        return S.One
    disc = (b * b - a * c).factor()
    if disc.is_positive:
        return S.One
    if disc.is_nonpositive and trace.is_positive:
        return S.NegativeOne
    return None


def _sqrt_complex(z: Expr) -> tuple[Expr, Expr]:
    """Computes the square root of a complex number.

    Unlike `sympy.sqrt`, it's correct for symbols that are potentially zero.
    See https://github.com/sympy/sympy/issues/29000 for more details.
    """
    x, y = z.as_real_imag()
    abs_z = sqrt(x * x + y * y)
    sgn = Piecewise((1, y >= 0), (-1, True))
    return (sqrt((abs_z + x) / 2), sgn * sqrt((abs_z - x) / 2))


def focal_axis_direction(conic: Conic) -> PointVector:
    """Returns the ideal point representing the direction of a conic's focal axis.

    The focal axis is also known as the transverse axis.

    Properties:
    - The focal axis is treated as an undirected line; its angle to the
      horizontal lies in the (-π/2, π/2] interval. For the full direction of a
      parabola, use [parabola_direction](#parabola.parabola_direction) instead.
    - Returns `[0, 0, 0]ᵀ` for circles and
      [circular conics](#conic_classes.is_circular), including zero-radius
      circles.
    - Point conics constructed by
      [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(ellipse)
      preserve the axis direction of the original real or imaginary ellipse.
    - Line pair conics constructed by
      [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(hyperbola)
      have no such property.
    - For crossing lines `l₁` and `l₂`, the focal axis of
      [line_pair_conic](#degenerate_conic.line_pair_conic)(l₁, l₂) and
      [angle_bisector](#line.angle_bisector)(l₁, l₂) point in the same
      direction.
    - The focal axis of parallel line pairs and double lines is parallel to
      the lines.

    *Formula*:
    [research/conic_properties/focus_directrix_eccentricity.py](../src/research/conic_properties/focus_directrix_eccentricity.py)

    *Related research* (trigonometric identities):
    [research/conic_properties/focal_axis_angle_identities.py](../src/research/conic_properties/focal_axis_angle_identities.py)
    """
    a, b, c, *_ = conic_coefficients(conic)
    norm_sign = ConicNormFactor(conic)
    x, y = _sqrt_complex(norm_sign * (a - c + 2 * I * b))
    return ideal_point(x, y)


def conjugate_axis_direction(conic: Conic) -> PointVector:
    """Returns the ideal point representing the direction of a conic's conjugate axis.

    It's equivalent to the ideal point on the directrix.

    Properties:
    - The conjugate axis is treated as an undirected line; its angle to the
      horizontal lies in the (0, π] interval.
    - Returns `[0, 0, 0]ᵀ` for circles and
      [circular conics](#conic_classes.is_circular), including zero-radius
      circles.
    - Point conics constructed by
      [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(ellipse)
      preserve the axis direction of the original real or imaginary ellipse.
    - Line pair conics constructed by
      [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(hyperbola)
      have no such property.
    """
    a, b, c, *_ = conic_coefficients(conic)
    norm_sign = ConicNormFactor(conic)
    x, y = _sqrt_complex(norm_sign * (a - c + 2 * I * b))
    return ideal_point(-y, x)
