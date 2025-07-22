from sympy import Expr, Matrix, Piecewise, sqrt

from lib.matrix import ConicMatrix, MaxEigenvalue, MinEigenvalue
from lib.point import PointToXY


def ConicFromCenterAndPoints(center, p1, p2, p3) -> Matrix:
    """Computes the conic section with the given center and perimeter points.

    May return
     - an ellipse;
     - a hyperbola;
     - a parallel line pair;
     - zero matrix if the solution is ambiguous, which happens when some of the
       (`center`, `pᵢ`, `pⱼ`) triples are collinear.
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
        ]
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


def ConicCenter(conic: Matrix):
    """Computes the center point of a conic.

    Formula: ChatGPT"""
    x, y, z = conic.row(0).cross(conic.row(1))
    return (x / z, y / z)


def SemiAxisLengths(conic: Matrix):
    """Computes the semi-axis lengths of a conic.

    Formula: ChatGPT
    """
    submatrix = conic[:2, :2]
    return (
        sqrt(-conic.det() / (MinEigenvalue(submatrix) * submatrix.det())),
        sqrt(-conic.det() / (MaxEigenvalue(submatrix) * submatrix.det())),
    )


def SemiMajorAxis(conic: Matrix) -> Expr:
    """Computes the semi-major axis length i.e. the center-vertex distance of
    a conic.

    It's infinity for parabolas, zero for degenerate conics, and imaginary for
    complex ellipses.
    """
    axes = SemiAxisLengths(conic)
    return Piecewise((axes[1], conic.det() > 0), (axes[0], True))


def SemiMinorAxis(conic: Matrix) -> Expr:
    """Computes the semi-minor axis length of a conic."""
    axes = SemiAxisLengths(conic)
    return Piecewise((axes[0], conic.det() > 0), (axes[1], True))
