from sympy import Matrix, sympify

ORIGIN = Matrix([0, 0, 1])


def IdealPoint(x, y):
    return Matrix([x, y, 0])


def PointToXY(point):
    assert len(point) in (2, 3)
    if len(point) == 2:
        return sympify(point)
    x, y, z = sympify(point)
    return (x / z, y / z)


def PointToVec3(point):
    assert len(point) in (2, 3)
    if len(point) == 2:
        return Matrix([point[0], point[1], 1])
    if point is Matrix:
        return point
    return Matrix(point)


def Centroid(*points):
    n = len(points)
    assert n > 0
    cx, cy = sympify((0, 0))
    for p in points:
        x, y = PointToXY(p)
        cx += x
        cy += y
    return (cx / n, cy / n)
