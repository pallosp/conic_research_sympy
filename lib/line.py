from sympy import Matrix

from lib.point import PointToVec3, PointToXY


def HorizontalLine(y):
    return Matrix([0, 1, -y])


def VerticalLine(x):
    return Matrix([-1, 0, x])


def LineBetween(point1, point2):
    return PointToVec3(point1).cross(PointToVec3(point2))


def ParallelLine(withLine, throughPoint):
    x, y = PointToXY(throughPoint)
    a, b, _ = withLine
    return Matrix([a, b, -a * x - b * y])


def AreParallel(line1, line2):
    """Tells whether line1 and line2 are parallel.

    Returns True if they are parallel, False if not, None if undecidable.
    Considers the ideal line parallel to everything.
    """
    a1, b1, _ = line1
    a2, b2, _ = line2
    return (a1 * b2 - a2 * b1).is_zero


IDEAL_LINE = Matrix([0, 0, 1])
X_AXIS = HorizontalLine(0)
Y_AXIS = VerticalLine(0)
