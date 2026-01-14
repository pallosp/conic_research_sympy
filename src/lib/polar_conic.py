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

from sympy import Expr, Matrix, atan2, cos, sin

from lib.circle import UNIT_CIRCLE
from lib.point import point_to_vec3

#: The circle at the origin with radius 1, in polar matrix form.
POLAR_UNIT_CIRCLE: Matrix = Matrix.eye(3)


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


def conic_from_polar_matrix(polar_conic: Matrix) -> Matrix:
    """Transforms a conic from polar to quadratic form.

    The algorithm is essentially applying the polar matrix as a projective
    transformation on the unit circle.
    """
    polar_adjugate = polar_conic.adjugate()
    return polar_adjugate.T * UNIT_CIRCLE * polar_adjugate
