from sympy import Matrix


def Circle(x, y, r_square):
    return Matrix(
        [
            [-1, 0, x],
            [0, -1, y],
            [x, y, r_square - x * x - y * y],
        ]
    )
