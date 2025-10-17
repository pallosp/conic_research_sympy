#!/usr/bin/env python3

"""Curvature sign computation for conics in polar form."""

from sympy import Determinant, Matrix, diff, pprint, sign, symbols

from lib.polar_conic import PointAtAngle

a, b, c, d, e, f, g, h, i, theta = symbols("a,b,c,d,e,f,g,h,i,theta")
polar_matrix = Matrix([[a, b, c], [d, e, f], [g, h, i]])
x, y, z = PointAtAngle(polar_matrix, theta)

# According to https://en.wikipedia.org/wiki/Curvature
# curvature = (x'y'' - y'x'') / (x'² + y'²)^(3/2)

# The conic has positive curvature at θ iff x'y'' - y'x'' > 0

xd = diff(x / z, theta)
xdd = diff(xd, theta)
yd = diff(y / z, theta)
ydd = diff(yd, theta)
discriminant = (xd * ydd - yd * xdd).factor().simplify()

# Multiplying with a full square doesn't change the sign but makes the formula
# work for points at infinity.
discriminant *= discriminant.as_numer_denom()[1].args[0] ** 4

# Un-evaluate the polar matrix determinant
discriminant = discriminant.subs(polar_matrix.det(), Determinant(polar_matrix))

print("\nSign of the polar conic's curvature at θ:\n")
pprint(sign(discriminant, evaluate=False))
