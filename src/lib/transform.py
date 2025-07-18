from sympy import cos, Matrix, sin


def TransformConic(conic: Matrix, transformation: Matrix) -> Matrix:
    return transformation.adjugate().T * conic * transformation.adjugate()


def TransformLine(line: Matrix, transformation: Matrix) -> Matrix:
    return transformation.adjugate().T * line


def Translate(dx, dy) -> Matrix:
    return Matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def Rotate(angle, x0=0, y0=0) -> Matrix:
    cos_angle = cos(angle)
    sin_angle = sin(angle)
    return Matrix(
        [
            [cos_angle, -sin_angle, x0 * (1 - cos_angle) + y0 * sin_angle],
            [sin_angle, cos_angle, y0 * (1 - cos_angle) - x0 * sin_angle],
            [0, 0, 1],
        ]
    )


def ScaleXY(scale_x, scale_y, x0=0, y0=0) -> Matrix:
    return Matrix(
        [
            [scale_x, 0, x0 * (1 - scale_x)],
            [0, scale_y, y0 * (1 - scale_y)],
            [0, 0, 1],
        ]
    )


def Scale(scale, x0=0, y0=0) -> Matrix:
    return ScaleXY(scale, scale, x0, y0)
