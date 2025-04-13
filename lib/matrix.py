from sympy import expand, Matrix, Not, Or, Piecewise, sqrt


def IsScalarMultiple(m1: Matrix, m2: Matrix) -> bool:
    """Tells whether two matrices are scalar multiples of each other."""
    if m1.shape != m2.shape:
        return False
    v1 = Matrix(list(m1))
    v2 = Matrix(list(m2))
    return expand(v1.dot(v1) * v2.dot(v2) - v1.dot(v2) ** 2) == 0


def MaxEigenvalue(symmetric_matrix2x2: Matrix):
    assert symmetric_matrix2x2.shape == (2, 2)
    a, b, b2, c = symmetric_matrix2x2
    assert b == b2
    return (a + c) / 2 + sqrt((a - c) ** 2 + 4 * b**2) / 2


def MinEigenvalue(symmetric_matrix2x2: Matrix):
    assert symmetric_matrix2x2.shape == (2, 2)
    a, b, b2, c = symmetric_matrix2x2
    assert b == b2
    return (a + c) / 2 - sqrt((a - c) ** 2 + 4 * b**2) / 2


def ConicMatrix(a, b, c, d, e, f):
    """3x3 symmetric matrix from the conic equation

            [a b d] [x]
    [x y 1] [b c e] [y] = 0
            [d e f] [1]

    Expanded form: ax² + 2bxy + cy² + 2dx + 2ey + f = 0
    """
    return Matrix([[a, b, d], [b, c, e], [d, e, f]])


def QuadraticForm(sym_matrix_3x3: Matrix, vector3: Matrix):
    """Quadratic form for a 3x3 symmetric matrix and a 3d column vector.

    Formula: vᵀ·M·v
    """
    return vector3.T * sym_matrix_3x3 * vector3


def SkewMatrix(vector3: Matrix):
    """Skew-symmetric matrix for a 3d vector."""
    x, y, z = vector3
    return Matrix(
        [
            [0, -z, y],
            [z, 0, -x],
            [-y, x, 0],
        ]
    )


def NonZeroRow(matrix: Matrix):
    """First non-zero row of the matrix or nan in case of zero matrix."""
    return Piecewise(
        *(
            (matrix.row(i), Or(*(Not(e.equals(0)) for e in matrix.row(i))))
            for i in range(matrix.rows)
        )
    )


def NonZeroCol(matrix: Matrix):
    """First non-zero column of the matrix or nan in case of zero matrix."""
    return Piecewise(
        *(
            (matrix.col(i), Or(*(Not(e.equals(0)) for e in matrix.col(i))))
            for i in range(matrix.cols)
        )
    )
