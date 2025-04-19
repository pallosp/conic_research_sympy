from sympy import Matrix

from lib.point import PointToVec3


def HorizontalLine(y):
    return Matrix([0, 1, -y])


def VerticalLine(x):
    return Matrix([-1, 0, x])


def LineBetween(point1, point2):
    return PointToVec3(point1).cross(PointToVec3(point2))


IDEAL_LINE = Matrix([0, 0, 1])
X_AXIS = HorizontalLine(0)
Y_AXIS = VerticalLine(0)
