#!/usr/bin/env python

from sympy import MatMul, factor, gcd, symbols

from lib.matrix import ConicMatrix
from lib.sympy_utils import EqChain
from research.util import println_indented

HORIZONTAL_LINE = "-" * 88

print(f"{HORIZONTAL_LINE}\n")

print("Conic matrix\n")

a, b, c, d, e, f = symbols("a b c d e f")
conic = ConicMatrix(a, b, c, d, e, f)
println_indented(conic)

print("Its adjugate\n")

adj_symbols = symbols("A_a B_a C_a D_a E_a F_a")
adj_matrix = ConicMatrix(*symbols("A_a B_a C_a D_a E_a F_a"))
adj_matrix_value = conic.adjugate()
println_indented(EqChain(symbols("adj"), adj_matrix, adj_matrix_value))

print("Adjugate of its adjugate\n")

double_adj = conic.adjugate().adjugate().applyfunc(factor)
common_factor = gcd(list(double_adj))
double_adj = MatMul(
    common_factor.subs(conic.det(), symbols("det")), (double_adj / common_factor),
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
