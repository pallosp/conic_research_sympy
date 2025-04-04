from sympy import Matrix


def ConicCenter(conic: Matrix):
    x, y, z = conic.row(0).cross(conic.row(1))
    return (x / z, y / z)
