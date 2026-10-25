#!/usr/bin/env python

from sympy import Matrix, Symbol, gcd, simplify, symbols

from lib.central_conic import conic_center
from lib.matrix import conic_matrix
from lib.transform import scale, transform_conic
from research.sympy_utils import println_indented

print("\nScaling a conic from its center by a factor of λ\n")

print("Original conic:\n")

a, b, c, d, e, f = symbols("a,b,c,d,e,f")
conic = conic_matrix(a, b, c, d, e, f)
println_indented(conic)

print("\nScaled conic - only the f coefficient differs:\n")

factor = symbols("lambda")
x, y = conic_center(conic).flat()
scaled_conic = transform_conic(conic, scale(factor, center=(x, y))).applyfunc(simplify)
scaled_conic = Matrix(scaled_conic / gcd(scaled_conic.flat()))
scaled_conic[8] = (scaled_conic[8] - f).factor() + f
scaled_conic[8] = scaled_conic[8].subs(conic.det(), Symbol("det"))

println_indented(scaled_conic)
