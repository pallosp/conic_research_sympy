from sympy import Matrix, nan, Piecewise, sqrt, Tuple


from lib.matrix import NonZeroCross, SkewMatrix


def LineXLine(line1: Matrix, line2: Matrix) -> Matrix:
    """Computes the intersection of two lines.

    Returns an ideal point if the lines are parallel, or a zero vector if they
    coincide.
    """
    return line1.cross(line2)


def ConicXLine(conic: Matrix, line: Matrix):
    """Intersects a conic with a line. Returns two points or NaN.

    Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 11.3
    """
    skew_matrix = SkewMatrix(line)
    m = skew_matrix.T * conic * skew_matrix
    a, b, c = line
    alpha = Piecewise(
        (sqrt(m[5] * m[7] - m[4] * m[8]) / a, a != 0),
        (sqrt(m[2] * m[6] - m[0] * m[8]) / b, b != 0),
        (sqrt(m[1] * m[3] - m[0] * m[4]) / c, c != 0),
    )
    if alpha == nan:
        return nan
    intersections = m + alpha * skew_matrix
    points = NonZeroCross(intersections)
    if isinstance(points, Tuple):
        return (points[0], points[1].T)
    return points
