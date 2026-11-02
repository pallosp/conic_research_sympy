#!/usr/bin/env python

"""Shows that the adjugate of the line pair conic `l₁ l₂ᵀ + l₂ l₁ᵀ` is `-x xᵀ`,
where `x = l₁ × l₂` is the intersection point of the two lines.
"""

from typing import cast

from sympy import Matrix, MatrixSymbol, Symbol, factor, symbols

from lib.degenerate_conic import line_pair_conic
from research.sympy_utils import eq_chain, println_indented

line1 = Matrix(symbols("a:c"))
line2 = Matrix(symbols("d:f"))

print("\nLine pair conic from lines l₁=(a,b,c) and l₂=(d,e,f):\n")

line_pair = line_pair_conic(line1, line2) * 2
println_indented(line_pair)

print("Its adjugate, factorized:\n")

adj = line_pair.adjugate().applyfunc(factor)
println_indented(adj)

print("Intersection of the lines:\n")

x = MatrixSymbol("x", 3, 1)
center = line1.cross(line2)
println_indented(eq_chain(x, center))
println_indented(eq_chain(x * x.T, center * center.T, -Symbol("adj(line pair)")))

assert cast("Matrix", adj).equals(-center * center.T)
