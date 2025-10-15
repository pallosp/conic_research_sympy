from collections.abc import Sequence

from sympy import Expr, Function, I, Matrix, Poly, Symbol, abc, sqrt

from lib.conic_classification import ConicNormFactor
from lib.matrix import NonZeroCross, QuadraticForm, SkewMatrix
from lib.point import PointToVec3, PointToXY


def ConicFromPoly(
    poly: Expr | Poly,
    *,
    x: Symbol = abc.x,
    y: Symbol = abc.y,
) -> Matrix:
    """Constructs the 3×3 symmetric matrix representation of a conic section
    from a two-variable quadratic polynomial in the form of
    ```
    ax² + bxy + cy² + dx + ey + f
    ```

    The resulting matrix is:
    ```
    [a, b/2, d/2]
    [b/2, c, e/2]
    [d/2, e/2, f]
    ```
    """
    poly = Poly(poly, x, y)
    a = poly.coeff_monomial(x * x)
    b = poly.coeff_monomial(x * y) / 2
    c = poly.coeff_monomial(y * y)
    d = poly.coeff_monomial(x) / 2
    e = poly.coeff_monomial(y) / 2
    f = poly.coeff_monomial(1)
    return Matrix([[a, b, d], [b, c, e], [d, e, f]])


def ConicThroughPoints(
    p1: Matrix | Sequence[Expr],
    p2: Matrix | Sequence[Expr],
    p3: Matrix | Sequence[Expr],
    p4: Matrix | Sequence[Expr],
    p5: Matrix | Sequence[Expr],
) -> Matrix:
    """Computes the conic that goes through the given points.

    Returns the conic matrix, or a zero matrix if the result is ambiguous.
    The result is unique if no 2 points coincide and no 4 points are collinear.
    The result is non-degenerate if no 3 points are collinear.

    *Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 10.1
    """
    p1, p2, p3, p4, p5 = [PointToVec3(p) for p in [p1, p2, p3, p4, p5]]
    g1 = p1.cross(p3)
    g2 = p2.cross(p4)
    h1 = p1.cross(p4)
    h2 = p2.cross(p3)
    g = g1 * g2.T + g2 * g1.T
    h = h1 * h2.T + h2 * h1.T
    return g * p5.dot(h1) * p5.dot(h2) - h * p5.dot(g1) * p5.dot(g2)


def ConicFromFocusAndDirectrix(
    focus: Matrix | Sequence[Expr],
    directrix: Matrix,
    eccentricity: Expr,
) -> Matrix:
    """Constructs a conic from its focus, directrix and eccentricity.

    *Formula*:
    [research/conic_from_focus_and_directrix.py](../src/research/conic_from_focus_and_directrix.py)
    """
    fx, fy = PointToXY(focus)
    m1 = directrix * directrix.T
    m2 = Matrix([[-1, 0, fx], [0, -1, fy], [fx, fy, -(fx**2) - fy**2]])
    return m1 * eccentricity**2 + m2 * (directrix[0] ** 2 + directrix[1] ** 2)


def Eccentricity(conic: Matrix) -> Expr:
    """Computes the eccentricity of a conic section.

    The result is
     - 0 for circles and imaginary circles;
     - (0..1) for other real ellipses;
     - imaginary for other imaginary ellipses;
     - 1 for parabolas;
     - &gt;1 for hyperbolas, in particular √2 for rectangular hyperbolas.

    In case of non-degenerate central conics, the eccentricity equals to the
    ratio of the center-focus distance
    ([LinearEccentricity](#central_conic.LinearEccentricity))
    and the center-vertex distance
    ([PrimaryRadius](#central_conic.PrimaryRadius)).

    The eccentricity of finite point conics constructed by shrinking an ellipse
    to zero size equals to that of the original ellipse.

    Crossing line pair conics have two different (generalized) focal axes, with
    two different corresponding eccentricity values. Evaluate
    `(Eccentricity(conic), Eccentricity(-conic))` to get both.

    *Formula*:
    https://en.wikipedia.org/wiki/Conic_section#Eccentricity_in_terms_of_coefficients<br>
    *Own research*:
    [research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)
    """
    a, b, c = conic[0], conic[1], conic[4]
    s = sqrt(((a - c) ** 2 + 4 * b**2).factor())
    norm_sign = ConicNormFactor(conic)
    return sqrt(2 * s / (s - norm_sign * (a + c)))


def FocalAxisDirection(conic: Matrix) -> Matrix:
    """Returns the direction of the focal axis a.k.a. transverse axis of the
    conic section, namely the ideal point in that direction.

    Properties:
    - Returns `[0, 0, 0]ᵀ` for circles and
      [circular conics](#conic_classification.IsCircular).
    - Point conics constructed by
      [ShrinkConicToZero](#central_conic.ShrinkConicToZero)(ellipse)
      preserve the axis direction of the original real or imaginary ellipse.
    - Line pair conics constructed by
      [ShrinkConicToZero](#central_conic.ShrinkConicToZero)(hyperbola)
      have no such property.
    - The focal axis of `LinePair(l1, l2)`, and `AngleBisector(l1, l2)` point to
      the same direction.

    *Formula*:
    [research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)
    """
    a, b, c = conic[0], conic[3], conic[4]
    norm_sign = ConicNormFactor(conic)
    x, y = sqrt(norm_sign * ((a - c) / 2 + b * I)).simplify().as_real_imag()
    return Matrix([x, y, 0])


class IdealPoints(Function):
    """Computes the ideal points on a conic section.

    Returns two points. Special cases:
     - For parabolas these are the same point.
     - For ellipses these are the complex conjugates of each other.
     - For symbolic conics returns an unevaluated `sympy.Function`.
    """

    @classmethod
    def eval(cls, conic: Matrix) -> tuple[Matrix, Matrix] | None:
        """Internal implementation. Call `IdealPoints(conic)` directly."""
        a, b, c = conic[0], conic[1], conic[4]
        disc = sqrt(b * b - a * c)
        cross = NonZeroCross(Matrix([[c, -b - disc, 0], [-b + disc, a, 0], [0, 0, 0]]))
        if isinstance(cross, NonZeroCross):
            return None
        return (cross[0], cross[1].T)


def ProjectiveConicCenter(conic: Matrix) -> Matrix:
    """Computes the generalized projective center of a conic.

    It's equivalent to [ConicCenter](#central_conic.ConicCenter) (returns a
    finite point) for
     - real and imaginary ellipses
     - hyperbolas
     - conics consisting of a single finite (Euclidean) point
     - crossing finite line pairs

    For parabolas returns the ideal point on it
    ([proof](../src/research/parabola_center.py))

    For other line pair conics and ideal point conics returns `[0, 0, 0]ᵀ`.
    """
    return conic.col(0).cross(conic.col(1))


def PolePoint(conic: Matrix, polar_line: Matrix) -> Matrix:
    """Computes the pole point of a conic with respect to the given polar line.

    If the conic is degenerate, i.e. it factors into `l₁` and `l₂` real or
    complex conjugate lines, the pole is
     - `[0, 0, 0]ᵀ` if `l₁`, `l₂`, and `polar_line` are concurrent;
     - the intersection of `l₁` and `l₂` otherwise.

    *Pole / polar identity*: `conic * pole_point = polar_line`<br>
    *Source*: https://en.wikipedia.org/wiki/Pole_and_polar#Calculating_the_pole_of_a_line
    """
    return conic.adjugate() * polar_line


def PolarLine(conic: Matrix, pole_point: Matrix | Sequence[Expr]) -> Matrix:
    """Computes the polar line of a conic with respect to the given pole point.

    *Pole / polar identity*: `conic * pole_point = polar_line`
    *Source*: https://en.wikipedia.org/wiki/Pole_and_polar#Calculating_the_pole_of_a_line
    """
    return conic * PointToVec3(pole_point)


def ConicContainsPoint(conic: Matrix, point: Matrix | Sequence[Expr]) -> bool | None:
    """Checks if a point lies on a conic.

    Returns None if undecidable.
    """
    return QuadraticForm(conic, PointToVec3(point)).expand().is_zero


def ConicContainsLine(conic: Matrix, line: Matrix) -> bool | None:
    """Checks if a line lies on a conic.

    Returns None if undecidable.

    *Formula*:
    [research/conic_line_containment.py](../src/research/conic_line_containment.py)
    """
    skew = SkewMatrix(line)
    return (skew * conic * skew).is_zero_matrix
