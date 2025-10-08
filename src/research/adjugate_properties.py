#!/usr/bin/env python

from sympy import Ge, MatMul, factor, gcd, symbols

from lib.matrix import ConicMatrix
from lib.sympy_utils import EqChain
from research.util import println_indented

HORIZONTAL_LINE = "-" * 80

print(f"{HORIZONTAL_LINE}\n")

print("Conic matrix\n")

a, b, c, d, e, f = symbols("a b c d e f")
conic = ConicMatrix(a, b, c, d, e, f)
println_indented(conic)

print("Its adjugate\n")

a_a, b_a, c_a, d_a, e_a, f_a = adj_symbols = symbols("A_a B_a C_a D_a E_a F_a")
adj_matrix = ConicMatrix(*symbols("A_a B_a C_a D_a E_a F_a"))
adj_matrix_value = conic.adjugate()
println_indented(EqChain(symbols("adj"), adj_matrix, adj_matrix_value))

print("Adjugate of its adjugate\n")

double_adj = conic.adjugate().adjugate().applyfunc(factor)
common_factor = gcd(list(double_adj))
det_value = conic.det()
det = symbols("det")
double_adj = MatMul(
    common_factor.subs(det_value, det),
    double_adj / common_factor,
)
println_indented(EqChain(symbols("double_adj"), adj_matrix.adjugate(), double_adj))

print("Special case: parabola (a⋅c - b² = Fₐ = 0)\n")

println_indented(
    EqChain(
        symbols("double_adj"),
        adj_matrix.adjugate().subs(adj_symbols[-1], 0),
        double_adj,
    ),
)

print(f"{HORIZONTAL_LINE}\n")

print("When the matrix is singular, the adjugate's diagonal elements are either all ≥0")
print("or all ≤0.\n")

println_indented(EqChain(adj_matrix.adjugate()[2, 2], det * f, 0))
println_indented(Ge(a_a * c_a, 0))
println_indented(Ge(conic.adjugate()[0] * conic.adjugate()[4], 0))

print("Similarly:\n")

println_indented(Ge(conic.adjugate()[0] * conic.adjugate()[8], 0))
println_indented(Ge(conic.adjugate()[4] * conic.adjugate()[8], 0))

print(f"{HORIZONTAL_LINE}\n")
