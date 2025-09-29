#!/usr/bin/env python

from sympy import Eq, MatMul, symbols

from lib.conic import ProjectiveConicCenter
from lib.matrix import ConicMatrix, QuadraticForm
from lib.sympy_utils import EqChain
from research.util import println_indented

HORIZONTAL_LINE = "-" * 100

a, b, c, d, e, f = symbols("a b c d e f")
parabola = ConicMatrix(a, b, c, d, e, f)

print(f"{HORIZONTAL_LINE}\n")

print("Conic matrix adjugate\n")

adj_symbols = ConicMatrix(*symbols("A_adj B_adj C_adj D_adj E_adj F_adj"))
adj_elements = parabola.adjugate()
aa, _, _, ab, ac, _, ad, ae, af = adj_elements
println_indented(EqChain(symbols("adj"), adj_elements, adj_symbols))

println_indented("The last row contains the generalized center point's coordinates.")

print(f"{HORIZONTAL_LINE}\n")

print("Lemma: parabola's projective center = parabola's ideal point\n")

center = ProjectiveConicCenter(parabola)
cx, cy, cz = center
println_indented(Eq(symbols("center"), center, evaluate=False))

print("It's an ideal point, because a·c - b² = 0\n")

print("It's on the parabola, because\n")

center = center.subs(a * c - b * b, 0)
center_quadratic_form = QuadraticForm(parabola, center).factor()
println_indented(EqChain(MatMul(center.T, parabola, center), center_quadratic_form, 0))

print(f"{HORIZONTAL_LINE}\n")

print("Parabola matrix determinant:\n")

det = symbols("det")
det_value = parabola.det().collect(f)
println_indented(Eq(det, det_value))

det_value = det_value.subs(a * c - b * b, 0)
println_indented(Eq(det, det_value))

print(f"{HORIZONTAL_LINE}\n")

print("Square sum of the center point's coordinates:\n")

center_square_sum = symbols("center.x") ** 2 + symbols("center.y") ** 2
center_square_sum_value = ad * ad + ae * ae
println_indented(Eq(center_square_sum, center_square_sum_value))

print("After expanding the formula and replacing b² with a·c:\n")

center_square_sum_value = center_square_sum_value.expand().subs(b * b, a * c).factor()
println_indented(
    EqChain(
        center_square_sum,
        center_square_sum_value,
        center_square_sum_value.replace(-det_value, -symbols("det")),
    ),
)

print("Corollary: sign(a + c) = -sign(det)\n")

print(f"{HORIZONTAL_LINE}\n")

print("d⋅det and e⋅det can also be expressed with the adjugate matrix elements.\n")

d_times_det = (d * parabola.det()).expand()
println_indented(Eq(d * det, d_times_det))

print("The c⋅d³ term hints at an (a⋅f - d²)⋅(b⋅e - c⋅d) subexpression.")
print("The remainder turns out to be the product of two other adj. matrix elements.\n")

d_times_det_factored = (d_times_det + ac * ad).factor() - ac * ad
println_indented(
    EqChain(
        d * det,
        d_times_det_factored,
        d_times_det_factored.subs(zip(adj_elements, adj_symbols, strict=True)),
    ),
)

print("Similarly\n")

e_times_det = (e * parabola.det()).expand()
e_times_det_factored = (e_times_det + aa * ae).factor() - aa * ae
println_indented(
    EqChain(
        e * det,
        e_times_det_factored,
        e_times_det_factored.subs(zip(adj_elements, adj_symbols, strict=True)),
    ),
)

print(f"{HORIZONTAL_LINE}\n")
