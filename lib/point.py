from sympy import Matrix

ORIGIN = Matrix([0, 0, 1])


def PointToXY(point):
    assert len(point) in (2, 3)
    return (
        (point[0], point[1])
        if len(point) == 2
        else (point[0] / point[2], point[1] / point[2])
    )
