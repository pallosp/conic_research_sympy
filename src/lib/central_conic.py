from collections.abc import Sequence

from sympy import Abs, Expr, Matrix, sqrt

from lib.conic import ConicNormFactor
from lib.matrix import ConicMatrix, MaxEigenvalue, MinEigenvalue
from lib.point import PointToXY


def ConicFromFociAndRadius(
    focus1: Matrix | Sequence[Expr],
    focus2: Matrix | Sequence[Expr],
    radius: Expr,
) -> Matrix:
    """Computes the ellipse or hyperbola with the given focus points and
    radius, i.e. center-vertex distance.

    If `radius` is negative, takes its absolute value.

    *Formula*:
    [research/conic_from_foci_and_radius.py](../src/research/conic_from_foci_and_radius.py)
    """
    fx1, fy1 = PointToXY(focus1)
    fx2, fy2 = PointToXY(focus2)

    # center
    cx, cy = (fx1 + fx2) / 2, (fy1 + fy2) / 2

    # center -> focus vector
    dx, dy = (fx2 - fx1) / 2, (fy2 - fy1) / 2

    a = dx**2 - radius**2
    b = dx * dy
    c = dy**2 - radius**2
    d = -cx * a - cy * b
    e = -cx * b - cy * c
    f = -cx * d - cy * e + (a * c - b * b)

    return ConicMatrix(a, b, c, d, e, f)


def ConicFromCenterAndPoints(
    center: Matrix | Sequence[Expr],
    p1: Matrix | Sequence[Expr],
    p2: Matrix | Sequence[Expr],
    p3: Matrix | Sequence[Expr],
) -> Matrix:
    """Computes the conic section with the given center and perimeter points.

    May return
     - an ellipse;
     - a hyperbola;
     - a parallel line pair;
     - zero matrix if the solution is ambiguous, which happens when some of the
       (`center`, `pᵢ`, `pⱼ`) triples are collinear.

    *Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)
    """
    x, y = PointToXY(center)
    x1, y1 = PointToXY(p1)
    x2, y2 = PointToXY(p2)
    x3, y3 = PointToXY(p3)
    x1, x2, x3 = x1 - x, x2 - x, x3 - x
    y1, y2, y3 = y1 - y, y2 - y, y3 - y
    m = Matrix(
        [
            [x1 * x1, x1 * y1, y1 * y1],
            [x2 * x2, x2 * y2, y2 * y2],
            [x3 * x3, x3 * y3, y3 * y3],
        ],
    )
    m1, m2, m3 = m.copy(), m.copy(), m.copy()
    m1[:, 0] = Matrix([1, 1, 1])
    m2[:, 1] = Matrix([1, 1, 1])
    m3[:, 2] = Matrix([1, 1, 1])
    a = m1.det()
    b = m2.det() / 2
    c = m3.det()
    d = -a * x - b * y
    e = -b * x - c * y
    f = -d * x - e * y - m.det()
    return ConicMatrix(a, b, c, d, e, f)


def ConicCenter(conic: Matrix) -> tuple[Expr, Expr]:
    """Computes the center point of a conic.

    *Formula*: [research/conic_center.py](../src/research/conic_center.py)
    """
    x, y, z = conic.col(0).cross(conic.col(1))
    return (x / z, y / z)


def SemiAxisLengths(conic: Matrix) -> tuple[Expr, Expr]:
    """Computes the semi-axis lengths of a conic in no specific order.

    To get the semi-focal or semi-transverse axis length (semi-major /
    semi-minor in case of ellipses), call
    [PrincipalRadius](#central_conic.PrincipalRadius) or
    [SecondaryRadius](#central_conic.SecondaryRadius), respectively.

    *Formula*: [research/conic_radii.py](../src/research/conic_radii.py)
    """
    submatrix = conic[:2, :2]
    return (
        sqrt(-conic.det() / (MinEigenvalue(submatrix) * submatrix.det())),
        sqrt(-conic.det() / (MaxEigenvalue(submatrix) * submatrix.det())),
    )


def PrimaryRadius(conic: Matrix) -> Expr:
    """Computes the center-vertex distance of a conic.

    This corresponds to the semi-major axis length of real ellipses. In case of
    complex ellipses however the focal axis is the shorter one in terms of
    absolute value.

    The returned value is:
     - a positive number for ellipses and hyperbolas;
     - infinity for parabolas;
     - an imaginary number for complex ellipses;
     - `nan` for ideal point conics;
     - 0 for the other degenerate conics.
    """
    a, _, _, b, c, _, _, _, _ = conic
    norm_sign = ConicNormFactor(conic)
    eigenvalue = (a + c + norm_sign * sqrt((a - c) ** 2 + 4 * b**2)) / 2
    return sqrt(-conic.det() / (eigenvalue * (a * c - b * b)))


def SecondaryRadius(conic: Matrix) -> Expr:
    """Computes the semi-conjugate axis length of a conic.

    This corresponds to the semi-minor axis length of real ellipses. In case of
    of complex ellipses however the conjugate axis is the longer one in terms
    absolute value. Hyperbolas intersect their conjugate axis at complex points,
    therefore the secondary radius will be a complex number.

    The returned value is:
     - a positive number for ellipses;
     - infinity for parabolas;
     - an imaginary number for hyperbolas and complex ellipses;
     - `nan` for ideal point conics;
     - 0 for the other degenerate conics.
    """
    a, _, _, b, c, _, _, _, _ = conic
    norm_sign = ConicNormFactor(conic)
    eigenvalue = (a + c - norm_sign * sqrt((a - c) ** 2 + 4 * b**2)) / 2
    return sqrt(-conic.det() / (eigenvalue * (a * c - b * b)))


def LinearEccentricity(conic: Matrix) -> Expr:
    """Computes the linear eccentricity of a conic section.

    The linear eccentricity is the distance between the center and a focus
    point. Special cases:
     - zero for circles and complex radius circles;
     - a real number for complex ellipses, since they still have real center
       and foci;
     - infinity for parabolas;
     - `nan` for parallel and coincident line pairs;
     - `nan` for conics containing the ideal line;
     - zero for all other degenerate conics.

    *Formula*: √|r₁²-r₂²| where r₁ and r₂ are the semi-axis lengths.
    """
    a, _, _, b, c, _, _, _, _ = conic
    l1, l2 = conic[:2, :2].eigenvals(multiple=True)
    return sqrt(Abs(conic.det() * (l1 - l2))) / Abs(a * c - b * b)


def ShrinkConicToZero(conic: Matrix) -> Matrix:
    """Scales a conic section from its center with a factor of zero.

    Turns hyperbolas to line pair conics consisting of their asymptotes, and
    ellipses to point conics.

    This transformation is only meaningful for central conics: for other
    conic types the result matrix will have infinite or `nan` elements.

    *Formula*:
    [research/scale_conic_from_center.py](../src/research/scale_conic_from_center.py
    """
    a, _, _, b, c, _, _, _, _ = conic
    return conic - Matrix.diag([0, 0, conic.det() / (a * c - b * b)])
