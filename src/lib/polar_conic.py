"""Utilities for conic sections in polar parametric form.

Polar conic sections are represented as parametric curves with angle θ, where
each point is the image of a point on the unit circle (`x² + y² = 1`) under a
projective transformation:

```
       [a b c]   [cos θ]
C(θ) = [d e f] * [sin θ]
       [g h i]   [  1  ]
```
"""

from sympy import cos, Matrix, sin

from lib.circle import UNIT_CIRCLE


POLAR_UNIT_CIRCLE = Matrix.eye(3)


def PointAtAngle(polar: Matrix, theta) -> Matrix:
    return polar * Matrix([cos(theta), sin(theta), 1])


def ConicFromPolarMatrix(polar: Matrix) -> Matrix:
    """Transforms a conic from polar to quadratic form.

    Algorithm: TransformConic(UNIT_CIRCLE, polar)
    """
    polar_adjugate = polar.adjugate()
    return polar_adjugate.T * UNIT_CIRCLE * polar_adjugate
