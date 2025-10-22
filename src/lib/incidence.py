from sympy import Matrix

from lib.point import PointToVec3


def AreCollinear(*points: Matrix) -> bool | None:
    """Tells whether n points are collinear.

    Returns `None` if undecidable.

    Algorithm:
     - n=3: https://en.wikipedia.org/wiki/Incidence_(geometry)#Collinearity
     - n>3: https://en.wikipedia.org/wiki/Gram_matrix
    """
    if len(points) <= 2:
        return True
    points_as_matrix = Matrix.hstack(*(PointToVec3(p) for p in points))
    if len(points) == 3:
        return points_as_matrix.det().is_zero
    return (points_as_matrix * points_as_matrix.T).det().is_zero
