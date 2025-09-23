#!/usr/bin/env python

from sympy import Eq, MatMul, symbols

from lib.conic import ProjectiveConicCenter
from lib.matrix import ConicMatrix, QuadraticForm
from research.util import println_indented

a, b, c, d, e, f = symbols("a b c d e f")
parabola = ConicMatrix(a, b, c, d, e, f)

print("\nLemma: parabola's projective center = parabola's ideal point\n")

center = ProjectiveConicCenter(parabola)
cx, cy, cz = center
println_indented(Eq(symbols("center"), center, evaluate=False))

print("It's an ideal point, because a·c - b² = 0\n")

print("It's on the parabola, because\n")

center = center.subs(a * c - b * b, 0)
center_quadratic_form = QuadraticForm(parabola, center).factor()
println_indented(
    Eq(MatMul(center.T, parabola, center), center_quadratic_form, evaluate=False),
)
