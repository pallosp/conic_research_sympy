#!/usr/bin/env python

from sympy import Abs, Symbol, sqrt, symbols

from lib.matrix import conic_matrix
from research.sympy_utils import eq_chain, println_indented

a, b, c, d, e, f = symbols("a b c d e f", real=True)
circle = conic_matrix(a, b, c, d, e, f)
det = Symbol("det", real=True)  # circle.det()
r = Symbol("r")

print("\n1. Radius from the conic matrix (det = determinant of the matrix):\n")

# See research/construction/director_circle.py
director_circle_radius = sqrt(-det * (a + c)) / Abs(a * c - b * b)
radius = director_circle_radius / sqrt(2)

println_indented(eq_chain(r, radius))

print("2. Circles: b=0 and c=a, keeping det symbolic:\n")

circle_radius = radius.subs({b: 0, c: a}).factor()

println_indented(eq_chain(r, circle_radius))

print("Formulae 1 and 2 both give 0 for finite non-circular point conics")
print("(det=0, a,c≠0, ac>b²).\n")

print("3. Circles only: expanding det after b=0 and c=a:\n")

circle_radius_expanded = (
    circle_radius.subs(det, circle.det())
    .subs({b: 0, c: a})
    .factor()
    .subs(Abs(a) / a**2, 1 / Abs(a))
)

println_indented(eq_chain(r, circle_radius_expanded))

print("Unlike formulae 1 and 2, this one doesn't hold for non-circular point conics.")
print("E.g. (x-1)² + 2(y-2)² = 0 has radius 0, but the formula gives 2√2.\n")
