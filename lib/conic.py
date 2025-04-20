from sympy import Matrix, Piecewise, sqrt

from lib.point import PointToVec3, PointToXY


def ConicThroughPoints(p1, p2, p3, p4, p5) -> Matrix:
    """Computes the conic that goes through the given points.

    The result is unique if no 2 points coincide and no 4 points are collinear.
    The result is non-degenerate if no 3 points are collinear.

    Returns the conic matrix, or a zero matrix if the result is ambiguous.

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
    focus: Matrix, directrix: Matrix, eccentricity
) -> Matrix:
    """Source: conic_from_focus_and_directrix.py"""
    fx, fy = PointToXY(focus)
    m1 = directrix * directrix.T
    m2 = Matrix([[-1, 0, fx], [0, -1, fy], [fx, fy, -(fx**2) - fy**2]])
    return m1 * eccentricity**2 + m2 * (directrix[0] ** 2 + directrix[1] ** 2)


def Eccentricity(conic: Matrix):
    """Eccentricity of the conic section.

    The result is ambiguous in case of degenerate conics: evaluate
    (Eccentricity(conic), Eccentricity(-conic)) to get both values.

    Source: https://en.wikipedia.org/wiki/Conic_section#Eccentricity_in_terms_of_coefficients
    """
    a, b, c = conic[0], conic[1], conic[4]
    s = sqrt(((a - c) ** 2 + 4 * b**2).factor())
    det_sign = Piecewise((1, conic.det() >= 0), (-1, True))
    return sqrt(2 * s / (s - det_sign * (a + c)))
