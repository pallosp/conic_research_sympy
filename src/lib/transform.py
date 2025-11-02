from collections.abc import Sequence

from sympy import Expr, Matrix, cos, sin

from lib.point import point_to_vec3


def transform_point(
    point: Matrix | Sequence[Expr],
    transformation: Matrix,
) -> Matrix:
    """Applies a projective transformation to a projective point.

    Returns a 3D column vector even if the input point is specified with only
    two coordinates.
    """
    return transformation * point_to_vec3(point)


def transform_line(line: Matrix, transformation: Matrix) -> Matrix:
    """Applies a projective transformation to a projective line."""
    return transformation.adjugate().T * line


def transform_conic(conic: Matrix, transformation: Matrix) -> Matrix:
    """Applies a projective transformation to a conic."""
    return transformation.adjugate().T * conic * transformation.adjugate()


def translate(dx: Expr, dy: Expr) -> Matrix:
    """Computes the transformation matrix for a 2D translation."""
    return Matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def rotate(angle: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
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


def reflect_to_line(axis: Matrix) -> Matrix:
    """Computes the transformation matrix for a reflection to a line.

    Returns a `nan` matrix if `axis` is the ideal line.

    *Formula*:
    [research/reflection_matrix.py](../src/research/reflection_matrix.py)
    """
    a, b, c = axis
    return Matrix(
        [
            [b**2 - a**2, -2 * a * b, -2 * a * c],
            [-2 * a * b, (a**2 - b**2), -2 * b * c],
            [0, 0, a**2 + b**2],
        ],
    ).applyfunc(lambda el: el / (a**2 + b**2))


def scale_xy(scale_x: Expr, scale_y: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
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


def scale(scale: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix:
    """Computes the projective transformation matrix for a uniform scaling
    transformation.
    """
    return scale_xy(scale, scale, x0, y0)


def _transformation_from_square_to_quad(target_points: Sequence[Matrix]) -> Matrix:
    """Computes the transformation that maps a square to a quadrilateral.

    It maps in particular
    - (1, 1) to `target_points[0]`
    - (-1, 1) to `target_points[1]`
    - (-1, -1) to `target_points[2]`
    - (1, -1) to `target_points[3]`

    where each target point is a 3d column vector with homogeneous coordinates.
    Returns a 3x3 projective transformation matrix.
    """
    x0, y0, z0 = target_points[0]
    x1, y1, z1 = target_points[1]
    x2, y2, z2 = target_points[2]
    x3, y3, z3 = target_points[3]
    coeffs = [
        x0 * x2 * (-y1 * z3 + y3 * z1)
        + x0 * x3 * (y1 * z2 - y2 * z1)
        + x1 * x2 * (y0 * z3 - y3 * z0)
        + x1 * x3 * (-y0 * z2 + y2 * z0),
        x0 * x1 * (y2 * z3 - y3 * z2)
        + x0 * x2 * (-y1 * z3 + y3 * z1)
        + x1 * x3 * (y0 * z2 - y2 * z0)
        + x2 * x3 * (-y0 * z1 + y1 * z0),
        x0 * x1 * (y2 * z3 - y3 * z2)
        + x0 * x3 * (y1 * z2 - y2 * z1)
        + x1 * x2 * (-y0 * z3 + y3 * z0)
        + x2 * x3 * (y0 * z1 - y1 * z0),
        y0 * y2 * (x1 * z3 - x3 * z1)
        + y0 * y3 * (-x1 * z2 + x2 * z1)
        + y1 * y2 * (-x0 * z3 + x3 * z0)
        + y1 * y3 * (x0 * z2 - x2 * z0),
        y0 * y1 * (-x2 * z3 + x3 * z2)
        + y0 * y2 * (x1 * z3 - x3 * z1)
        + y1 * y3 * (-x0 * z2 + x2 * z0)
        + y2 * y3 * (x0 * z1 - x1 * z0),
        y0 * y1 * (-x2 * z3 + x3 * z2)
        + y0 * y3 * (-x1 * z2 + x2 * z1)
        + y1 * y2 * (x0 * z3 - x3 * z0)
        + y2 * y3 * (-x0 * z1 + x1 * z0),
        z0 * z2 * (-x1 * y3 + x3 * y1)
        + z0 * z3 * (x1 * y2 - x2 * y1)
        + z1 * z2 * (x0 * y3 - x3 * y0)
        + z1 * z3 * (-x0 * y2 + x2 * y0),
        z0 * z1 * (x2 * y3 - x3 * y2)
        + z0 * z2 * (-x1 * y3 + x3 * y1)
        + z1 * z3 * (x0 * y2 - x2 * y0)
        + z2 * z3 * (-x0 * y1 + x1 * y0),
        z0 * z1 * (x2 * y3 - x3 * y2)
        + z0 * z3 * (x1 * y2 - x2 * y1)
        + z1 * z2 * (-x0 * y3 + x3 * y0)
        + z2 * z3 * (x0 * y1 - x1 * y0),
    ]
    return Matrix(3, 3, coeffs)


def transformation_from_samples(
    source_points: Sequence[Matrix | Sequence[Expr]],
    target_points: Sequence[Matrix | Sequence[Expr]],
) -> Matrix:
    """Computes the transformation that maps one quadrilateral to another.

    Takes 4 projective source and 4 projective target points. Both the source
    and the target quadrilaterals must be non-degenerate. Returns a 3x3
    projective transformation matrix.

    *Algorithm*:
    [research/homography_from_samples.py](../src/research/homography_from_samples.py)
    """
    if len(source_points) != 4 or len(target_points) != 4:
        raise ValueError("Exactly 4 source and 4 target points are required")

    source_points = [point_to_vec3(p) for p in source_points]
    target_points = [point_to_vec3(p) for p in target_points]

    t1 = _transformation_from_square_to_quad(source_points)
    t2 = _transformation_from_square_to_quad(target_points)
    return t1.inv() * t2
