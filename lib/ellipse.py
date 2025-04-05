from sympy import cos, Matrix, sin


def Ellipse(center_x, center_y, r1, r2, r1_angle=0):
    """Source: ellipse_from_params.py"""
    a = -((r1 * sin(r1_angle)) ** 2) - (r2 * cos(r1_angle)) ** 2
    b = (r1**2 - r2**2) * sin(r1_angle) * cos(r1_angle)
    c = -((r1 * cos(r1_angle)) ** 2) - (r2 * sin(r1_angle)) ** 2
    d = -a * center_x - b * center_y
    e = -b * center_x - c * center_y
    f = r1**2 * r2**2 - d * center_x - e * center_y
    return Matrix([[a, b, d], [b, c, e], [d, e, f]])
