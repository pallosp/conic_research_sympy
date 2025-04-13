from sympy import Matrix, Piecewise, sqrt

from lib.matrix import NonZeroCol, NonZeroRow, SkewMatrix


def ConicXLine(conic: Matrix, line: Matrix):
    """Source: Jürgen Richter-Gebert, Projective Geometry, section 11.3"""
    skew_matrix = SkewMatrix(line)
    m = skew_matrix.T * conic * skew_matrix
    a, b, c = line
    alpha = Piecewise(
        (sqrt(m[5] * m[7] - m[4] * m[8]) / a, a != 0),
        (sqrt(m[2] * m[6] - m[0] * m[8]) / b, b != 0),
        (sqrt(m[1] * m[3] - m[0] * m[4]) / c, c != 0),
    )
    intersections = m + alpha * skew_matrix
    return (NonZeroCol(intersections), NonZeroRow(intersections).T)
