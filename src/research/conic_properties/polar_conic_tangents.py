#!/usr/bin/env python3

from sympy import Matrix, cos, expand, gcd, sin, symbols

from lib.conic import polar_line
from lib.polar_conic import conic_from_polar_matrix, point_at_angle
from research.sympy_utils import println_indented

polar_matrix = Matrix(3, 3, symbols("a,b,c,d,e,f,g,h,i"))
angle = symbols("alpha")
point = point_at_angle(polar_matrix, angle)
conic = conic_from_polar_matrix(polar_matrix)

tangent = polar_line(conic, point)
tangent = tangent.applyfunc(lambda el: el.factor().collect([cos(angle), sin(angle)]))
tangent /= gcd(list(tangent))

print("\nTangent to a conic in polar representation:\n")

println_indented(tangent)

print("which equals to adjugate(P)ᵀ·[-cos α, -sin α, 1]ᵀ.\n")

expected_tangent = polar_matrix.adjugate().T * Matrix([-cos(angle), -sin(angle), 1])
assert expand(tangent) == expand(expected_tangent)

print("  OK\n")
