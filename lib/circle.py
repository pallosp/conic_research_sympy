from lib.matrix import ConicMatrix
from lib.point import ORIGIN, PointToXY


def Circle(center, r):
    x, y = PointToXY(center)
    return ConicMatrix(-1, 0, -1, x, y, r * r - x * x - y * y)


UNIT_CIRCLE = Circle(ORIGIN, 1)
