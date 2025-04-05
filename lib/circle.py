from sympy import Matrix


def Circle(x, y, r):
    return Matrix(
        [
            [-1, 0, x],
            [0, -1, y],
            [x, y, r * r - x * x - y * y],
        ]
    )


def UnitCircle():
    return Circle(0, 0, 1)
