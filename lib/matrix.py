from sympy import Function, Matrix, nan, sqrt


def IsNonZeroMultiple(m1: Matrix | list, m2: Matrix | list) -> bool:
    """Tells whether two matrices are non-zero scalar multiples of each other.

    Treats lists as column vectors.
    """
    if not isinstance(m1, Matrix):
        m1 = Matrix(m1)
    if not isinstance(m2, Matrix):
        m2 = Matrix(m2)
    if m1.shape != m2.shape:
        return False
    if m1.is_zero_matrix != m2.is_zero_matrix:
        return False
    v1 = m1.reshape(len(m1), 1)
    v2 = m2.reshape(len(m2), 1)
    return (v1.dot(v1) * v2.dot(v2) - v1.dot(v2) ** 2).equals(0)


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


def IsDefinite(symmetric_matrix: Matrix) -> bool:
    """Checks if a symmetric matrix is positive or negative definite."""
    assert symmetric_matrix.is_symmetric()

    e0 = symmetric_matrix[0]
    positive = None
    if e0.is_positive:
        positive = True
    elif e0.is_zero:
        return False
    elif e0.is_negative:
        positive = False
    else:
        return None

    for i in range(2, symmetric_matrix.rows + 1):
        minor_det = symmetric_matrix[:i, :i].det()
        must_be_positive = i % 2 == 0 or positive
        if must_be_positive:
            if minor_det.is_positive:
                continue
            if minor_det.is_nonpositive:
                return False
            return None
        else:
            if minor_det.is_negative:
                continue
            if minor_det.is_nonnegative:
                return False
            return None

    return True
