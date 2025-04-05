from sympy import Matrix, Piecewise, sqrt

from lib.matrix import MaxEigenvalue, MinEigenvalue


def ConicCenter(conic: Matrix):
    """Source: ChatGPT"""
    x, y, z = conic.row(0).cross(conic.row(1))
    return (x / z, y / z)


def AxisLengths(conic: Matrix):
    """Source: ChatGPT"""
    submatrix = conic[:2, :2]
    return (
        sqrt(-conic.det() / (MinEigenvalue(submatrix) * submatrix.det())),
        sqrt(-conic.det() / (MaxEigenvalue(submatrix) * submatrix.det())),
    )


def SemiMajorAxis(conic: Matrix):
    axes = AxisLengths(conic)
    return Piecewise((axes[1], conic.det() > 0), (axes[0], True))


def SemiMinorAxis(conic: Matrix):
    axes = AxisLengths(conic)
    return Piecewise((axes[0], conic.det() > 0), (axes[1], True))
