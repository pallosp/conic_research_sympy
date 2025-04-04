from sympy import Matrix


def TransformConic(conic, matrix):
    return matrix.adjugate().T * conic * matrix.adjugate()


def Translate(dx, dy):
    return Matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])
