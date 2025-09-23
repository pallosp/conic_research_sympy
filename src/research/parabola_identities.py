#!/usr/bin/env python

from sympy import Eq, Lt, MatMul, symbols

from lib.conic import ProjectiveConicCenter
from lib.matrix import ConicMatrix, QuadraticForm
from research.util import println_indented

HORIZONTAL_LINE = "-" * 88

a, b, c, d, e, f = symbols("a b c d e f")
parabola = ConicMatrix(a, b, c, d, e, f)

print(f"\n{HORIZONTAL_LINE}\n")

print("Lemma: parabola's projective center = parabola's ideal point\n")

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

print(f"{HORIZONTAL_LINE}\n")

print("Determinant:\n")

det_value = parabola.det().collect(f)
println_indented(Eq(symbols("det"), det_value))

det_value = det_value.subs(a * c - b * b, 0)
println_indented(Eq(symbols("det"), det_value))

print(f"{HORIZONTAL_LINE}\n")

print("Sign of a + c:\n")

center_square_sum = symbols("center.x") ** 2 + symbols("center.y") ** 2
center_square_sum_value = cx * cx + cy * cy
println_indented(Eq(center_square_sum, center_square_sum_value))

center_square_sum_value = center_square_sum_value.expand()
println_indented(Eq(center_square_sum, center_square_sum_value))

print("After replacing b² with a·c:\n")

center_square_sum_value = center_square_sum_value.subs(b * b, a * c).factor()
println_indented(Eq(center_square_sum, center_square_sum_value))

center_square_sum_value = center_square_sum_value.replace(-det_value, -symbols("det"))
println_indented(Lt(0, center_square_sum_value))

print(f"{HORIZONTAL_LINE}\n")
