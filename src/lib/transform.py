from collections.abc import Sequence

from sympy import Expr, Matrix, cos, sin, sympify


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
        ],
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
        ],
    )


def Scale(scale: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
    """Computes the projective transformation matrix for a uniform scaling
    transformation.
    """
    return ScaleXY(scale, scale, x0, y0)


def TransformationFromSamples(
    source_points: Sequence[tuple[Expr, Expr]],
    target_points: Sequence[tuple[Expr, Expr]],
) -> Matrix:
    """Computes the transformation matrix that maps one quadrilateral to
    another.

    *Algorithm*:
    https://franklinta.com/2014/09/08/computing-css-matrix3d-transforms/
    """
    if len(source_points) != 4 or len(target_points) != 4:
        raise ValueError("Exactly 4 source and 4 target points are required")
    (x0, y0), (x1, y1), (x2, y2), (x3, y3) = source_points
    (u0, v0), (u1, v1), (u2, v2), (u3, v3) = target_points
    m = Matrix(
        [
            [x0, y0, 1, 0, 0, 0, -u0 * x0, -u0 * y0],
            [0, 0, 0, x0, y0, 1, -v0 * x0, -v0 * y0],
            [x1, y1, 1, 0, 0, 0, -u1 * x1, -u1 * y1],
            [0, 0, 0, x1, y1, 1, -v1 * x1, -v1 * y1],
            [x2, y2, 1, 0, 0, 0, -u2 * x2, -u2 * y2],
            [0, 0, 0, x2, y2, 1, -v2 * x2, -v2 * y2],
            [x3, y3, 1, 0, 0, 0, -u3 * x3, -u3 * y3],
            [0, 0, 0, x3, y3, 1, -v3 * x3, -v3 * y3],
        ],
    )
    uv = Matrix([u0, v0, u1, v1, u2, v2, u3, v3])
    coefficients = [*list(m.inv() * uv), sympify(1)]
    return Matrix(3, 3, coefficients)
