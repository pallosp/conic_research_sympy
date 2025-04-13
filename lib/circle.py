from lib.matrix import ConicMatrix


def Circle(x, y, r):
    return ConicMatrix(-1, 0, -1, x, y, r * r - x * x - y * y)


def UnitCircle():
    return Circle(0, 0, 1)
