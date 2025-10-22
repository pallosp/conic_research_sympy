#!/usr/bin/env python3

from sympy import MatMul, Matrix, factor, gcd, pprint, symbols

from lib.polar_conic import conic_from_polar_matrix

a, b, c, d, e, f, g, h, i = symbols("a,b,c,d,e,f,g,h,i")
polar_matrix = Matrix([[a, b, c], [d, e, f], [g, h, i]])

conic = conic_from_polar_matrix(polar_matrix)
center = conic.row(0).cross(conic.row(1)).T
center /= -gcd(center[0], center[1])
center = center.applyfunc(factor)

assert center == polar_matrix * Matrix([g, h, -i])

print("\nCenter point of a polar conic:\n")
pprint(MatMul(polar_matrix, Matrix([g, h, -i])))
