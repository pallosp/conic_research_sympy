from sympy import sqrt

from lib.point import PointToXY, PointToXYZ


def PointPointDistance(point1, point2):
    """Signed distance between two points.

    Infinity if one of them is an ideal point, NaN if both."""
    x1, y1, z1 = PointToXYZ(point1)
    x2, y2, z2 = PointToXYZ(point2)
    return sqrt((x2 * z1 - x1 * z2) ** 2 + (y2 * z1 - y1 * z2) ** 2) / (z1 * z2)


def PointLineDistance(point, line):
    """Signed distance between a point and a line.

    Infinity for an Euclidean point and an ideal line.
    NaN for ideal points.
    """
    x, y = PointToXY(point)
    a, b, c = line
    return (a * x + b * y + c) / sqrt(a * a + b * b)
