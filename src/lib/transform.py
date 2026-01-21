from collections.abc import Sequence

from sympy import Expr, Matrix, cos, sin

from lib.point import ORIGIN, point_to_vec3, point_to_xy


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


def transform_polar_conic(polar_conic: Matrix, transformation: Matrix) -> Matrix:
    """Applies a projective transformation to a conic in polar representation."""
    return transformation * polar_conic


def translate(by: Matrix | Sequence[Expr]) -> Matrix:
    """Computes the transformation matrix for a 2D translation."""
    dx, dy = by
    return Matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])


def rotate(angle_radians: Expr, around: Matrix | Sequence[Expr] = ORIGIN) -> Matrix:
    """Computes the transformation matrix for a rotation around a point."""
    x0, y0 = point_to_xy(around)
    cos_angle = cos(angle_radians)
    sin_angle = sin(angle_radians)

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
    [research/transformation/reflection_matrix.py](../src/research/transformation/reflection_matrix.py)
    """
    a, b, c = axis
    return Matrix(
        [
            [b**2 - a**2, -2 * a * b, -2 * a * c],
            [-2 * a * b, (a**2 - b**2), -2 * b * c],
            [0, 0, a**2 + b**2],
        ],
    ).applyfunc(lambda el: el / (a**2 + b**2))


def scale_xy(
    scale_x: Expr,
    scale_y: Expr,
    center: Sequence[Expr] | Matrix = ORIGIN,
) -> Matrix:
    """Computes the projective transformation matrix for scaling along the x-
    and y-axes.
    """
    x0, y0 = point_to_xy(center)
    return Matrix(
        [
            [scale_x, 0, x0 * (1 - scale_x)],
            [0, scale_y, y0 * (1 - scale_y)],
            [0, 0, 1],
        ],
    )


def scale(scale: Expr, center: Sequence[Expr] | Matrix = ORIGIN) -> Matrix:
    """Computes the projective transformation matrix for a uniform scaling
    transformation.
    """
    return scale_xy(scale, scale, center=center)


def homography_from_samples(
    source_points: Sequence[Matrix | Sequence[Expr]],
    target_points: Sequence[Matrix | Sequence[Expr]],
) -> Matrix:
    """Computes the transformation that maps one quadrilateral to another.

    Takes 4 projective source and 4 projective target points. Both the source
    and the target quadrilaterals must be non-degenerate. Returns a 3x3
    projective transformation matrix.

    *Research*:
    [research/transformation/homography_from_samples.py](../src/research/transformation/homography_from_samples.py)
    """
    if len(source_points) != 4 or len(target_points) != 4:
        raise ValueError("Exactly 4 source and 4 target points are required")

    # 3x4 matrices of the source and target points
    s = Matrix.hstack(*[point_to_vec3(p) for p in source_points])
    t = Matrix.hstack(*[point_to_vec3(p) for p in target_points])

    # Determinants of s and t without columns 0, 1 or 2
    s0, s1, s2 = [s[:, [j for j in range(4) if j != i]].det() for i in range(3)]
    t0, t1, t2 = [t[:, [j for j in range(4) if j != i]].det() for i in range(3)]

    # Step 1: transform the source points to (1,0,0), (0,1,0), (0,0,1), (1,1,1).
    # Step 2: transform (1,0,0), (0,1,0), (0,0,1), (1,1,1) to the target points.
    #
    # The intermediade points could be anything, but this specific selection
    # results in very simple transformation matrices for both steps.
    return t[:, :3] * Matrix.diag([t0 / s0, t1 / s1, t2 / s2]) * s[:, :3].inv()
