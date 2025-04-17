from sympy import cos, sin

from lib.matrix import ConicMatrix
from lib.point import PointToXY


def Ellipse(center, r1, r2, *, r1_angle=None, r1_direction=None):
    """Source: ellipse_from_params.py"""
    assert r1_angle is None or r1_direction is None
    center_x, center_y = PointToXY(center)
    axis_dir_x, axis_dir_y = (
        r1_direction or (1, 0) if r1_angle is None else (cos(r1_angle), sin(r1_angle))
    )
    a = -((r1 * axis_dir_y) ** 2) - (r2 * axis_dir_x) ** 2
    b = (r1**2 - r2**2) * axis_dir_x * axis_dir_y
    c = -((r1 * axis_dir_x) ** 2) - (r2 * axis_dir_y) ** 2
    d = -a * center_x - b * center_y
    e = -b * center_x - c * center_y
    f = (r1**2 * r2**2) * (axis_dir_x**2 + axis_dir_y**2) - d * center_x - e * center_y
    return ConicMatrix(a, b, c, d, e, f)
