from sympy import sqrt

from lib.point import PointToXY


def PointLineDistance(point, line):
    """Signed distance between a point and a line."""
    x, y = PointToXY(point)
    a, b, c = line
    return (a * x + b * y + c) / sqrt(a * a + b * b)
