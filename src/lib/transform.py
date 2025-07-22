from sympy import Expr, cos, Matrix, sin


def TransformConic(conic: Matrix, transformation: Matrix) -> Matrix:
    """Applies a projective transformation on a conic."""
    return transformation.adjugate().T * conic * transformation.adjugate()


def TransformLine(line: Matrix, transformation: Matrix) -> Matrix:
    """Applies a projective transformation on a projective line."""
    return transformation.adjugate().T * line


def Translate(dx: Expr, dy: Expr) -> Matrix:
    """Computes the transformation matrix for a 2d translation."""
    return Matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def Rotate(angle: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
    """Computes the transformation matrix for a rotation around a point."""
    cos_angle = cos(angle)
    sin_angle = sin(angle)
    return Matrix(
        [
            [cos_angle, -sin_angle, x0 * (1 - cos_angle) + y0 * sin_angle],
            [sin_angle, cos_angle, y0 * (1 - cos_angle) - x0 * sin_angle],
            [0, 0, 1],
        ]
    )


def ScaleXY(scale_x: Expr, scale_y: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
    """Computes the projective transformation matrix for scaling along the x-
    and y-axes.
    """
    return Matrix(
        [
            [scale_x, 0, x0 * (1 - scale_x)],
            [0, scale_y, y0 * (1 - scale_y)],
            [0, 0, 1],
        ]
    )


def Scale(scale: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
    """Computes the projective transformation matrix for a uniform scaling
    transformation.
    """
    return ScaleXY(scale, scale, x0, y0)
