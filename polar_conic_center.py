#!/usr/bin/env python3

from sympy import MatMul, Matrix, gcd, pprint, symbols

from lib.polar_conic import ConicFromPolarMatrix

a, b, c, d, e, f, g, h, i = symbols("a,b,c,d,e,f,g,h,i")
polar_matrix = Matrix([[a, b, c], [d, e, f], [g, h, i]])

conic = ConicFromPolarMatrix(polar_matrix)
center = conic.row(0).cross(conic.row(1)).T
center /= -gcd(center[0], center[1])
for index in range(3):
    center[index] = center[index].factor()

assert center == polar_matrix * Matrix([g, h, -i])

pprint("\nCenter point of a polar conic:\n\n")
pprint(MatMul(polar_matrix, Matrix([g, h, -i])))
