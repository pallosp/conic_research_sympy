from typing import cast

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
from sympy.core.logic import fuzzy_and, fuzzy_or
from sympy.core.parameters import global_parameters

from lib.matrix import conic_coefficients
from lib.point import ideal_point
from lib.type_aliases import Conic, PointVector, Simplifier


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
    - *Double ideal line*: The factor `-sign(f)` makes the conic equation
      `f·w²` negative at every finite point, like for finite double lines.
    - *Crossing lines, other conics containing the ideal line*: No preferred
      normalization exists. The factor is 1.

    The table lists the resulting definiteness of the normalized matrix and of
    its upper left 2×2 block. It also lists the formulas that give the factor.
    In the formulas, `det` is the determinant of `C`. `a`, `c` and `f` are the
    diagonal elements of `C`. Circles count as ellipses.

    | conic               | matrix        | 2×2 block       | factor                   |
    |---------------------|---------------|-----------------|--------------------------|
    | ellipse             | indefinite    | neg. definite   | `sign(det) = -sign(a+c)` |
    | imaginary ellipse   | pos. definite | pos. definite   | `sign(det) = sign(a+c)`  |
    | parabola            | indefinite    | neg. semidef. ¹ | `sign(det) = -sign(a+c)` |
    | hyperbola           | indefinite    | indefinite      | `sign(det)`              |
    | finite point conic  | neg. semidef. | neg. definite   | `-sign(a+c)`             |
    | ideal point conic   | neg. semidef. | neg. semidef.   | `-sign(a+c)`             |
    | parallel lines      | indefinite    | neg. semidef.   | `-sign(a+c)`             |
    | finite double line  | neg. semidef. | neg. semidef.   | `-sign(a+c)`             |
    | double ideal line   | neg. semidef. | zero            | `-sign(f)`               |
    | crossing lines      | indefinite    | indefinite      | 1                        |
    | finite ∪ ideal line | indefinite    | zero            | 1                        |
    | zero matrix         | zero          | zero            | 1                        |

    ¹ Follows from being positive at the focus
    ([proof](../src/research/conic_properties/normalized_parabola.py)).

    `C.det() * ConicNormFactor(C) == Abs(C.det())` holds for all conic types.

    Takes an optional `simplifier` callback that simplifies the polynomials
    before their signs get determined. Symbolic conics may remain unevaluated
    if undecidable. Since the result doesn't remember the simplifier, use
    `expr.doit(simplifier=...)` to retry the evaluation later.
    """

    # Implies is_real = True, is_integer = True, and is_nonzero = True.
    is_odd = True

    # mypy expects the return type to be ConicNormFactor
    def __new__(  # type: ignore[misc]
        cls,
        conic: Conic,
        *,
        simplifier: Simplifier = lambda expr: expr,
        evaluate: bool | None = None,
    ) -> Expr:
        """Creates the factor, simplifying the polynomials with `simplifier`."""
        # SymPy's `eval` hook would only receive the Basic args, not
        # `simplifier`, so the evaluation happens here.
        if evaluate is None:
            evaluate = global_parameters.evaluate
        if not evaluate:
            return super().__new__(cls, conic, evaluate=False)

        det = simplifier(conic.det())
        if det.is_nonzero:  # implies det.is_real
            return sign(det)  # ±1 if det's sign is known

        if det.is_zero:
            factor = _degenerate_conic_norm_factor(conic, simplifier)
            if factor is not None:
                return factor

        # It's undecidable whether the conic is degenerate. The result is
        # valid if sign(det) and the degenerate conic's factor agree.
        if det.is_nonnegative and _degenerate_conic_norm_factor(conic, simplifier) == 1:
            return S.One

        if (
            det.is_nonpositive
            and _degenerate_conic_norm_factor(conic, simplifier) == -1
        ):
            return S.NegativeOne

        return super().__new__(cls, conic, evaluate=False)

    def doit(self, **hints: object) -> Expr:
        """Retries the evaluation, simplifying with the `simplifier` hint."""
        conic = cast("Conic", self.args[0])
        if hints.get("deep", True):
            conic = conic.doit(**hints)
        simplifier = cast("Simplifier", hints.get("simplifier", lambda expr: expr))
        return ConicNormFactor(conic, simplifier=simplifier)

    def _eval_Abs(self) -> Integer:  # noqa: N802 (hook name looked up by Abs)
        return S.One

    def _eval_power(self, exponent: Expr) -> Expr | None:
        if exponent.is_even:
            return S.One
        if exponent.is_odd:
            return self
        return None


def _degenerate_conic_norm_factor(
    conic: Conic,
    simplifier: Simplifier,
) -> Expr | None:
    """Computes `ConicNormFactor` for degenerate conics.

    The factor is -1 in two cases, and 1 otherwise:
    - The discriminant `b² - a·c` is ≤ 0 and the trace `a+c` is positive. This
      covers point conics, parallel lines and finite double lines. A nonzero
      semidefinite 2×2 block has the sign of its trace, so the factor makes
      the block negative semidefinite.
    - The conic is a double ideal line with `f > 0`.
    """
    a, b, c, d, e, f = conic_coefficients(conic)
    trace = a + c
    disc = simplifier(b * b - a * c)
    flip = fuzzy_or(
        [
            fuzzy_and([trace.is_positive, disc.is_nonpositive]),
            fuzzy_and([f.is_positive, *(x.is_zero for x in (a, b, c, d, e))]),
        ],
    )
    if flip is None:
        return None
    return S.NegativeOne if flip else S.One


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
