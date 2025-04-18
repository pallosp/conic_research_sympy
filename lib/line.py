from sympy import Matrix


def HorizontalLine(y):
    return Matrix([0, 1, -y])


def VerticalLine(x):
    return Matrix([-1, 0, x])


IDEAL_LINE = Matrix([0, 0, 1])
X_AXIS = HorizontalLine(0)
Y_AXIS = VerticalLine(0)
