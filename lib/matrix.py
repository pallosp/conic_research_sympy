from sympy import expand, Function, Matrix, nan, sqrt


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


def QuadraticForm(sym_matrix: Matrix, vector: Matrix):
    """Quadratic form for a nxn symmetric matrix and a n-element column vector.

    Formula: vᵀ·M·v
    """
    return vector.dot(vector.T * sym_matrix)


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


class NonZeroCross(Function):
    @classmethod
    def eval(cls, matrix: Matrix):
        all_zero = True
        for i in range(len(matrix)):
            el = matrix[i]
            if el.is_nonzero:
                return (matrix.col(i % matrix.cols), matrix.row(i // matrix.cols))
            if all_zero and not el.is_zero:
                all_zero = False
        if all_zero:
            return nan
