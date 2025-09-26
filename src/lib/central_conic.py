from collections.abc import Sequence

from sympy import Expr, Function, Matrix, sqrt

from lib.conic_classification import IsFiniteConic
from lib.matrix import ConicMatrix, MaxEigenvalue, MinEigenvalue
from lib.point import PointToXY


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
    """Computes the semi-axis lengths of a conic.

    *Formula*: [research/conic_radii.py](../src/research/conic_radii.py)
    """
    submatrix = conic[:2, :2]
    return (
        sqrt(-conic.det() / (MinEigenvalue(submatrix) * submatrix.det())),
        sqrt(-conic.det() / (MaxEigenvalue(submatrix) * submatrix.det())),
    )


class SemiMajorAxis(Function):
    """Computes the semi-major axis length i.e. the center-vertex distance of
    a conic.

    The returned value is:
     - a positive number for ellipses and hyperbolas;
     - infinity for parabolas;
     - an imaginary number for complex ellipses;
     - nan for ideal point conics;
     - zero for the other degenerate conics.

    Returns an unevaluated `sympy.Function` if we can't tell which axis is
    longer.
    """

    @classmethod
    def eval(cls, conic: Matrix) -> Expr | None:
        """Internal implementation. Call `SemiMajorAxis(conic)` directly."""
        axes = SemiAxisLengths(conic)
        finite = IsFiniteConic(conic)
        if finite is True:
            if conic[0].is_negative:
                return axes[1]
            if conic[0].is_nonnegative:
                return axes[0]
        if finite is False:
            det = conic.det()
            if det.is_nonnegative:
                return axes[1]
            if det.is_negative:
                return axes[0]
        return None


class SemiMinorAxis(Function):
    """Computes the semi-minor axis length of a conic.

    The returned value is:
     - a positive number for ellipses;
     - infinity for parabolas;
     - an imaginary number for hyperbolas and complex ellipses;
     - nan for ideal point conics;
     - zero for the other degenerate conics.

    Returns an unevaluated `sympy.Function` if we can't tell which axis is
    shorter.
    """

    @classmethod
    def eval(cls, conic: Matrix) -> Expr | None:
        """Internal implementation. Call `SemiMinorAxis(conic)` directly."""
        axes = SemiAxisLengths(conic)
        finite = IsFiniteConic(conic)
        if finite is True:
            if conic[0].is_negative:
                return axes[0]
            if conic[0].is_nonnegative:
                return axes[1]
        if finite is False:
            det = conic.det()
            if det.is_nonnegative:
                return axes[0]
            if det.is_negative:
                return axes[1]
        return None
