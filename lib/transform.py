from sympy import cos
from sympy import Matrix
from sympy import sin


def TransformConic(conic, matrix):
    return matrix.adjugate().T * conic * matrix.adjugate()


def Translate(dx, dy):
    return Matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def Rotate(angle, x0=0, y0=0):
    cos_angle = cos(angle)
    sin_angle = sin(angle)
    return Matrix(
        [
            [cos_angle, -sin_angle, x0 * (1 - cos_angle) + y0 * sin_angle],
            [sin_angle, cos_angle, y0 * (1 - cos_angle) - x0 * sin_angle],
            [0, 0, 1],
        ]
    )
