#!/usr/bin/env python3

from sympy import MatMul, Matrix, factor, gcd, symbols

from lib.conic import projective_conic_center
from lib.polar_conic import conic_from_polar_matrix
from research.sympy_utils import println_indented

polar_matrix = Matrix(3, 3, symbols("a,b,c,d,e,f,g,h,i"))

conic = conic_from_polar_matrix(polar_matrix)
center = projective_conic_center(conic)
center /= -gcd(center[0], center[1])
center = center.applyfunc(factor)

print("\nCenter point of a conic in polar representation:\n")

println_indented(MatMul(polar_matrix, (polar_matrix.inv() * center).applyfunc(factor)))
