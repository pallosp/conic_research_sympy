#!/usr/bin/env python

from typing import cast

from sympy import Expr, Ge, MatMul, factor, gcd, symbols

from lib.matrix import conic_matrix
from research.sympy_utils import eq_chain, println_indented

HORIZONTAL_LINE = "-" * 80

print(HORIZONTAL_LINE)

print("Conic matrix\n")

a, b, c, d, e, f = symbols("a b c d e f")
conic = conic_matrix(a, b, c, d, e, f)
println_indented(conic)

print("Its adjugate\n")

a_a, b_a, c_a, d_a, e_a, f_a = adj_symbols = symbols("A_a B_a C_a D_a E_a F_a")
adj_matrix = conic_matrix(*symbols("A_a B_a C_a D_a E_a F_a"))
adj_matrix_value = conic.adjugate()
println_indented(eq_chain(symbols("adj"), adj_matrix, adj_matrix_value))

print("Adjugate of its adjugate\n")

double_adj = conic.adjugate().adjugate().applyfunc(factor)
# gcd() is typed to also return polynomial domain elements; with expressions it
# returns an Expr (could be fixed upstream with an overload).
common_factor = cast("Expr", gcd(double_adj.flat()))
det_value = conic.det()
det = symbols("det")
double_adj_factored = MatMul(
    common_factor.subs(det_value, det),
    double_adj / common_factor,
)
println_indented(
    eq_chain(symbols("double_adj"), adj_matrix.adjugate(), double_adj_factored)
)

print("Special case: parabola (a⋅c - b² = Fₐ = 0)\n")

println_indented(
    eq_chain(
        symbols("double_adj"),
        adj_matrix.adjugate().subs(adj_symbols[-1], 0),
        double_adj_factored,
    ),
)

print(HORIZONTAL_LINE)

print("When the matrix is singular, the adjugate's diagonal elements are either all ≥0")
print("or all ≤0.\n")

println_indented(eq_chain(adj_matrix.adjugate()[2, 2], det * f, 0))
println_indented(Ge(a_a * c_a, 0))
println_indented(Ge(conic.adjugate()[0] * conic.adjugate()[4], 0))

print("Similarly:\n")

println_indented(Ge(conic.adjugate()[0] * conic.adjugate()[8], 0))
println_indented(Ge(conic.adjugate()[4] * conic.adjugate()[8], 0))

print(HORIZONTAL_LINE)
