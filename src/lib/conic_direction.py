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

    *Geometric meaning*: The equation of the normalized conic `N = ConicNormFactor(C)·C`
    is positive at the points [inside](#incidence.is_inside) the conic and negative at
    the points [outside](#incidence.is_outside). For real non-degenerate conics, the
    inside points are the ones on the focus side(s), through which no real tangent
    passes.

    *Linear algebra meaning*: The first applicable rule determines the factor.

    1. If `det(C) ≠ 0`, the factor `sign(det(C))` makes `det(N)` positive.
    2. If the upper left 2×2 block `[[a, b], [b, c]]` is semidefinite and nonzero
       (`b² ≤ a·c` and `a+c ≠ 0`), the factor `-sign(a+c)` makes it negative
       semidefinite in `N`. If `C` is semidefinite, `N` becomes negative
       semidefinite too.
    3. If `C` is semidefinite and nonzero, the factor `-sign(f)` where `f = C[2, 2]`
       makes `N` negative semidefinite. After rule 2, this only applies when `f` is the
       only nonzero element of `C`.
    4. Otherwise, the factor is 1.

    The table lists the definiteness of `N` and of its upper left 2×2 block.
    Circles count as ellipses.

    | conic               | N             | N[:2, :2]       | factor                   |
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

    ¹ Follows from the normalized conic equation being positive at the focus
    ([proof](../src/research/conic_properties/normalized_parabola.py)).

    Consequences of the rules above:

    1. Makes the inside points of non-degenerate conics evaluate to positive.<br>
       `det(C)·ConicNormFactor(C) = |det(C)|` holds for all conic types.
    2. Keeps the normalization continuous: as `ε → 0`, the normalized parabolas
       `y² = εx + k` converge to the normalized parallel lines (`k > 0`), double line
       (`k = 0`) or ideal point conic (`k < 0`).
    3. Normalizes the double ideal line the same way as the finite double lines
       (negative semidefinite).
    4. Crossing lines, finite ∪ ideal line and the zero matrix have no preferred
       normalization: `ConicNormFactor(C) = ConicNormFactor(-C) = 1`.

    `ConicNormFactor` takes an optional `simplifier` callback that simplifies the
    polynomials before their signs get determined. Symbolic conics may remain
    unevaluated if undecidable. Since the unevaluated function doesn't remember the
    simplifier, use `expr.doit(simplifier=...)` to retry the evaluation later.
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
