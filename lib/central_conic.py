from sympy import Matrix, sqrt

from lib.matrix import MaxEigenvalue, MinEigenvalue


def ConicCenter(conic: Matrix):
    """Source: ChatGPT"""
    x, y, z = conic.row(0).cross(conic.row(1))
    return (x / z, y / z)


def SemiMajorAxis(conic: Matrix):
    """Source: ChatGPT"""
    submatrix = conic[:2, :2]
    return sqrt(-conic.det() / (MinEigenvalue(submatrix) * submatrix.det()))


def SemiMinorAxis(conic: Matrix):
    """Source: ChatGPT"""
    submatrix = conic[:2, :2]
    return sqrt(-conic.det() / (MaxEigenvalue(submatrix) * submatrix.det()))
