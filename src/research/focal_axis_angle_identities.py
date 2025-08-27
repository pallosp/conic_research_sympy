#!/usr/bin/env python

"""Computes how to express the trigonometric functions of the focal axis angle
from a conic matrix.

The starting point is the formula in focus_directrix_eccentricity.py:

    θ = atan2(2 * b * sign(det), (a - c) * sign(det)) / 2
"""

from sympy import Abs, Eq, atan2, cos, pprint, sign, simplify, sin, symbols

a, b, c, theta = symbols("a,b,c,theta", real=True)
det = symbols("det", real=True, nonzero=True)

print("\nFocal axis direction:\n")

angle = atan2(2 * b * sign(det), (a - c) * sign(det)) / 2
pprint(Eq(theta, angle))
print()

pprint(Eq(cos(theta * 2), simplify(cos(angle * 2))).subs(det / Abs(det), sign(det)))
print()

pprint(Eq(sin(theta * 2), simplify(sin(angle * 2))).subs(det / Abs(det), sign(det)))
print()

# cos²(θ) = (1 + cos(2θ)) / 2
cos_angle_square = simplify((1 + cos(angle * 2)) / 2).subs(det / Abs(det), sign(det))
pprint(Eq(cos(theta) ** 2, cos_angle_square))
print()

# sin²(θ) = (1 - cos(2θ)) / 2
sin_angle_square = simplify((1 - cos(angle * 2)) / 2).subs(det / Abs(det), sign(det))
pprint(Eq(sin(theta) ** 2, sin_angle_square))
print()
