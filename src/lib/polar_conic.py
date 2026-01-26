"""Utilities for conic sections in polar parametric form.

Polar conic sections are represented as parametric curves with angle θ, where
each point is the image of a point on the unit circle (`x² + y² = 1`) under a
projective transformation:

```
       ⎡a  b  c⎤   ⎡cos θ⎤
C(θ) = ⎢d  e  f⎥ * ⎢sin θ⎥
       ⎣g  h  i⎦   ⎣  1  ⎦
```
"""

from collections.abc import Sequence
from enum import Enum

from sympy import Expr, I, Matrix, atan2, cos, sign, sin, sqrt

from lib.central_conic import conic_center, primary_radius, secondary_radius
from lib.circle import UNIT_CIRCLE
from lib.conic_direction import focal_axis_direction
from lib.point import point_to_vec3

#: The circle at the origin with radius 1, in polar matrix form.
POLAR_UNIT_CIRCLE: Matrix = Matrix.eye(3)


class PolarOrigin(Enum):
    """Specifies which point of a conic in polar form corresponds to angle 0."""

    #: A vertex of the conic.
    VERTEX = 0

    #: A co-vertex of the conic.
    COVERTEX = 1

    #: An ideal point of the conic.
    IDEAL_POINT = 2

    #: An endpoint of the horizontal diameter.
    HORIZONTAL = 3

    #: An endpoint of the vertical diameter.
    VERTICAL = 4


def point_at_angle(polar_conic: Matrix, theta: Expr) -> Matrix:
    """Computes the coordinates of the projective point on a polar conic
    corresponding to a certain angle.
    """
    return polar_conic * Matrix([cos(theta), sin(theta), 1])


def angle_at_point(polar_conic: Matrix, point: Matrix | Sequence[Expr]) -> Expr:
    """Computes the polar angle corresponding to a point on a polar conic.

    The result is unspecified if the point is not on the conic.
    """
    point = point_to_vec3(point)
    x_times_cos_a, x_times_sin_a, _x = polar_conic.adjugate() * point
    return atan2(x_times_sin_a, x_times_cos_a)


def tangent_at_angle(polar_conic: Matrix, angle_radians: Expr) -> Matrix:
    """Computes the tangent line to a polar conic at the given angle.

    *Formula*:
    [research/conic_properties/polar_conic_tangents.py](../src/research/conic_properties/polar_conic_tangents.py)
    """
    adj = polar_conic.adjugate()
    return adj.T * Matrix([cos(angle_radians), sin(angle_radians), -1])


def curvature_sign_at_angle(polar_conic: Matrix, angle_radians: Expr) -> Matrix:
    """Tells which direction a polar conic turns at an angle.

    - *Positive*: the curve turns left (counterclockwise).
    - *Negative*: the curve turns right (clockwise).
    - *Zero*: the curve point at the angle is an ideal point.

    *Formula*:
    [research/conic_properties/polar_conic_curvature_sign.py](../src/research/conic_properties/polar_conic_curvature_sign.py)
    """
    g, h, i = polar_conic.row(2)
    det = polar_conic.det()
    return sign((g * cos(angle_radians) + h * sin(angle_radians) + i) * det)


def conic_from_polar_matrix(polar_conic: Matrix) -> Matrix:
    """Transforms a conic from polar to quadratic form.

    The algorithm is essentially applying the polar matrix as a projective
    transformation on the unit circle.
    """
    polar_adjugate = polar_conic.adjugate()
    return polar_adjugate.T * UNIT_CIRCLE * polar_adjugate


def ellipse_to_polar_matrix(
    ellipse: Matrix,
    start: PolarOrigin = PolarOrigin.HORIZONTAL,
) -> Matrix:
    """Converts an ellipse to a polar conic matrix representation.

    The resulting polar form has the following properties:
     - The vertices and covertices are separated by an angular distance of π/2.
     - The secants between α and α+π go through the ellipse center.
     - The homogeneous z-coordinates of all points on the curve are equal to 1.
     - The curve goes counterclockwise. To reverse its direction, negate the
       middle column of the matrix.
     - The point corresponding to angle 0 is determined by the `start` parameter.

    *Formula*:
    [research/construction/polar_ellipse.py](../src/research/construction/polar_ellipse.py)
    """
    a, _, _, b, c, _, d, e, _ = ellipse
    disc = a * c - b * b
    cx = (b * e - c * d) / disc
    cy = (b * d - a * e) / disc

    if start == PolarOrigin.HORIZONTAL:
        det = ellipse.det()
        t = sqrt(-det / a)
        return Matrix(
            [
                [t / sqrt(disc), b * t * sign(det) / disc, cx],
                [0, -a * t * sign(det) / disc, cy],
                [0, 0, 1],
            ]
        )

    if start == PolarOrigin.VERTICAL:
        det = ellipse.det()
        t = sqrt(-det / c)
        return Matrix(
            [
                [0, c * t * sign(det) / disc, cx],
                [t / sqrt(disc), -b * t * sign(det) / disc, cy],
                [0, 0, 1],
            ]
        )

    if start == PolarOrigin.VERTEX:
        fd = focal_axis_direction(ellipse)
        cos_a, sin_a, _ = fd / fd.norm()
        r1 = primary_radius(ellipse)
        r2 = secondary_radius(ellipse)

        return Matrix(
            [
                [r1 * cos_a, -r2 * sin_a, cx],
                [r1 * sin_a, r2 * cos_a, cy],
                [0, 0, 1],
            ]
        )

    raise ValueError("Unsupported PolarOrigin value")


def hyperbola_to_polar_matrix(
    hyperbola: Matrix,
    start: PolarOrigin = PolarOrigin.VERTEX,
) -> Matrix:
    """Converts a hyperbola to a polar conic matrix representation.

    The resulting polar form has the following properties:
     - The vertices and ideal points are separated by an angular distance of π/2.
     - The secants between α and α+π go through the center point.
     - The point corresponding to angle 0 is determined by the `start` parameter.

    *Formula*:
    [research/construction/polar_hyperbola.py](../src/research/construction/polar_hyperbola.py)
    """
    if start == PolarOrigin.VERTEX:
        fd = focal_axis_direction(hyperbola)
        cos_a, sin_a, _ = fd / fd.norm()
        cx, cy = conic_center(hyperbola)
        r1 = primary_radius(hyperbola)
        r2 = secondary_radius(hyperbola)

        return Matrix(
            [
                [cx, -I * r2 * sin_a, r1 * cos_a],
                [cy, I * r2 * cos_a, r1 * sin_a],
                [1, 0, 0],
            ]
        )

    raise ValueError("Unsupported PolarOrigin value")
