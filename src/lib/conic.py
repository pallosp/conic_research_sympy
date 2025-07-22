from sympy import Expr, Function, I, Matrix, Piecewise, Poly, Symbol, Tuple, abc, sqrt

from lib.matrix import NonZeroCross, SkewMatrix
from lib.point import PointToVec3, PointToXY


def ConicFromPoly(
    poly: Expr | Poly,
    x: Symbol = abc.x,
    y: Symbol = abc.y,
) -> Matrix:
    """
    Constructs the 3×3 symmetric matrix representation of a conic section
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


def ConicThroughPoints(p1, p2, p3, p4, p5) -> Matrix:
    """Computes the conic that goes through the given points.

    Returns the conic matrix, or a zero matrix if the result is ambiguous.
    The result is unique if no 2 points coincide and no 4 points are collinear.
    The result is non-degenerate if no 3 points are collinear.

    Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 10.1
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
    focus: Matrix, directrix: Matrix, eccentricity: Expr
) -> Matrix:
    """Constructs a conic from its focus, directrix and eccentricity.

    Formula: `research/conic_from_focus_and_directrix.py`
    """
    fx, fy = PointToXY(focus)
    m1 = directrix * directrix.T
    m2 = Matrix([[-1, 0, fx], [0, -1, fy], [fx, fy, -(fx**2) - fy**2]])
    return m1 * eccentricity**2 + m2 * (directrix[0] ** 2 + directrix[1] ** 2)


def Eccentricity(conic: Matrix):
    """Eccentricity of the conic section.

    The result is ambiguous in case of degenerate conics: evaluate
    (Eccentricity(conic), Eccentricity(-conic)) to get both values.

    Source: https://en.wikipedia.org/wiki/Conic_section#Eccentricity_in_terms_of_coefficients
    Own research: focus_directrix_eccentricity.py
    """
    a, b, c = conic[0], conic[1], conic[4]
    s = sqrt(((a - c) ** 2 + 4 * b**2).factor())
    det_sign = Piecewise((1, conic.det() >= 0), (-1, True))
    return sqrt(2 * s / (s - det_sign * (a + c)))


def AxisDirection(conic: Matrix) -> Matrix:
    """Returns the direction of the major axis of the conic section,
    namely the ideal point in that direction.

    The result is ambiguous in case of degenerate conics: evaluate
    (AxisDirection(conic), AxisDirection(-conic)) to get both values.
    """
    a, b, c = conic[0], conic[3], conic[4]
    sign = Piecewise((1, conic.det() >= 0), (-1, True))
    x, y = sqrt(sign * (a - c) / 2 + sign * b * I).simplify().as_real_imag()
    return Matrix([x, y, 0])


class SplitToLines(Function):
    """Splits a degenerate conic into two lines.

    In case of point conics the lines will be complex conjugates.

    Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 11.1
    """

    @classmethod
    def eval(cls, conic: Matrix):
        adj = conic.adjugate()
        A, C, F = adj.diagonal()

        if A.is_nonzero:
            conic = conic + SkewMatrix(adj.col(0) / sqrt(-A))
        elif C.is_nonzero:
            conic = conic + SkewMatrix(adj.col(1) / sqrt(-C))
        elif F.is_nonzero:
            conic = conic + SkewMatrix(adj.col(2) / sqrt(-F))
        elif not (A.is_zero and C.is_zero and F.is_zero):
            return None

        cross = NonZeroCross(conic)
        if isinstance(cross, Tuple):
            return (cross[0], cross[1].T)
        return None


class IdealPoints(Function):
    """Computes the ideal points of a conic section.

    Always returns two points. For parabolas these are the same point.
    For ellipses these are the complex conjugates of each other.
    """

    @classmethod
    def eval(cls, conic: Matrix):
        a, b, c = conic[0], conic[1], conic[4]
        disc = sqrt(b * b - a * c)
        cross = NonZeroCross(Matrix([[c, -b - disc, 0], [-b + disc, a, 0], [0, 0, 0]]))
        if isinstance(cross, Tuple):
            return (cross[0], cross[1].T)
        return None
