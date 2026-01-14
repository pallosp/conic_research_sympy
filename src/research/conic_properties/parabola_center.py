#!/usr/bin/env python

from sympy import MatMul, symbols

from lib.conic import projective_conic_center
from lib.matrix import conic_matrix, quadratic_form
from research.sympy_utils import eq_chain, println_indented

print("\nLemma: parabola's projective center = parabola's ideal point\n")

a, b, c, d, e, f = symbols("a b c d e f")
parabola = conic_matrix(a, b, c, d, e, f)
center = projective_conic_center(parabola)

print("The center point is the conic matrix adjugate's last column:\n")

println_indented(eq_chain(symbols("adj"), parabola.adjugate()))
println_indented(eq_chain(symbols("center"), center))

print("It's not a zero vector, because then the adjugate matrix would be singular.")

print("It's an ideal point, because a·c - b² = 0.")

print("It's on the parabola, because\n")

center = center.subs(a * c - b * b, 0)
center_quadratic_form = quadratic_form(parabola, center).factor()
println_indented(eq_chain(MatMul(center.T, parabola, center), center_quadratic_form, 0))
