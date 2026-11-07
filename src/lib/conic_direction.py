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

    Multiplying `C` by this factor gives the normalized matrix. The normalized
    matrix has the following properties:

    - *Non-degenerate conics*: The conic equation is positive at the focus
      point(s), i.e. `[fx fy 1]·C·[fx fy 1]ᵀ > 0`.
    - *Finite point conics*: The normalized matrix is negative semidefinite,
      i.e. the conic equation is ≤0 at every point. The upper left 2×2
      block is a principal submatrix, so it's negative semidefinite too. Its
      discriminant is negative, so it's even negative definite. A negative
      definite matrix has a negative trace, which is `factor·(a+c)` here.
      Hence the factor is `-sign(a+c)`.
    - *Ideal point conics, parallel lines and double lines*: These conics have
      a zero discriminant, like parabolas. So their upper left 2×2 block is
      singular. A singular symmetric 2×2 matrix is semidefinite. The factor
      `-sign(a+c)` makes the block negative semidefinite, as it is for
      parabolas
      ([proof](../src/research/conic_properties/normalized_parabola.py)).<br>
      The choice of factor keeps the normalization continuous. As `ε → 0`, the
      normalized parabolas `y² = εx + k` converge to the normalized parallel
      lines (`k > 0`), double line (`k = 0`) or ideal point conic (`k < 0`).
      For ideal point conics, the whole normalized matrix becomes negative
      semidefinite too, as for finite point conics.
    - *Crossing lines, conics containing the ideal line*: No preferred
      normalization exists. The factor is 1.

    The table lists the resulting definiteness of the normalized matrix and of
    its upper left 2×2 block. It also lists the formulas that give the factor.
    In the formulas, `det` is the determinant of `C`. `a` and `c` are the first
    two diagonal elements of `C`. Circles count as ellipses.

    | conic              | matrix        | 2×2 block       | factor                   |
    |--------------------|---------------|-----------------|--------------------------|
    | ellipse            | indefinite    | neg. definite   | `sign(det) = -sign(a+c)` |
    | imaginary ellipse  | pos. definite | pos. definite   | `sign(det) = sign(a+c)`  |
    | parabola           | indefinite    | neg. semidef. ¹ | `sign(det) = -sign(a+c)` |
    | hyperbola          | indefinite    | indefinite      | `sign(det)`              |
    | finite point conic | neg. semidef. | neg. definite   | `-sign(a+c)`             |
    | ideal point conic  | neg. semidef. | neg. semidef.   | `-sign(a+c)`             |
    | parallel lines     | indefinite    | neg. semidef.   | `-sign(a+c)`             |
    | double line        | neg. semidef. | neg. semidef.   | `-sign(a+c)`             |
    | crossing lines     | indefinite    | indefinite      | 1                        |
    | other degenerate ² | —             | zero            | 1                        |

    ¹ Follows from being positive at the focus
    ([proof](../src/research/conic_properties/normalized_parabola.py)).<br>
    ² Conics containing the ideal line, including the zero matrix.

    `C.det() * ConicNormFactor(C) == Abs(C.det())` holds for all conic types.

    Symbolic conics may remain unevaluated if undecidable.
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
