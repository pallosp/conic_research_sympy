#!/usr/bin/env python

from sympy import Matrix, Symbol, gcd, simplify, symbols

from lib.central_conic import ConicCenter
from lib.matrix import ConicMatrix
from lib.transform import Scale, TransformConic
from research.util import println_indented

print("\nScaling a conic from its center by a factor of λ\n")

print("Original conic:\n")

a, b, c, d, e, f = symbols("a,b,c,d,e,f")
conic = ConicMatrix(a, b, c, d, e, f)
println_indented(conic)

print("\nScaled conic - only the f coefficient differs:\n")

factor = symbols("lambda")
x, y = ConicCenter(conic)
scaled_conic = simplify(TransformConic(conic, Scale(factor, x, y)))
scaled_conic = Matrix(scaled_conic / gcd(list(scaled_conic)))
scaled_conic[8] = (scaled_conic[8] - f).factor() + f
scaled_conic[8] = scaled_conic[8].subs(conic.det(), Symbol("det"))

println_indented(scaled_conic)
