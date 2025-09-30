#!/usr/bin/env python

from sympy import MatMul, symbols

from lib.conic import ProjectiveConicCenter
from lib.matrix import ConicMatrix, QuadraticForm
from lib.sympy_utils import EqChain
from research.util import println_indented

print("\nLemma: parabola's projective center = parabola's ideal point\n")

a, b, c, d, e, f = symbols("a b c d e f")
parabola = ConicMatrix(a, b, c, d, e, f)
center = ProjectiveConicCenter(parabola)

print("The center point is the conic matrix adjugate's last column:\n")

println_indented(EqChain(symbols("adj"), parabola.adjugate()))
println_indented(EqChain(symbols("center"), center))

print("It's not a zero vector, because then the adjugate matrix would be singular.")

print("It's an ideal point, because a·c - b² = 0.")

print("It's on the parabola, because\n")

center = center.subs(a * c - b * b, 0)
center_quadratic_form = QuadraticForm(parabola, center).factor()
println_indented(EqChain(MatMul(center.T, parabola, center), center_quadratic_form, 0))
