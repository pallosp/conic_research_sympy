from sympy import Matrix, Piecewise, sqrt

from lib.point import PointToXY


def ConicFromFocusAndDirectrix(focus: Matrix, directrix: Matrix, eccentricity):
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
