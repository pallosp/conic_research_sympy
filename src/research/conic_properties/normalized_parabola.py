#!/usr/bin/env python

"""Proves two lemmas about conic normalization.

- The upper left 2×2 block of a normalized parabola is negative semidefinite.
- In order to make the upper left 2×2 block of a conic with exactly one ideal
  point negative semidefinite, the normalization factor must be -sign(a + c).
"""

from typing import cast

from sympy import Abs, Le, MatMul, Matrix, MatrixSymbol, Symbol, sign, symbols

from lib.conic import conic_from_focus_and_directrix
from lib.matrix import quadratic_form
from research.equations import eq_chain
from research.pretty_print import HORIZONTAL_LINE, println_indented

print(
    "\nLemma: the upper left 2×2 block of a normalized parabola's matrix is\n"
    "negative semidefinite.\n"
)

print(
    "Every parabola is a non-zero multiple of the one constructed from its\n"
    "focus F and directrix [a, b, c]ᵀ with eccentricity 1:\n"
)

a, b, c, fx, fy = symbols("a b c fx fy", real=True)
k = Symbol("k", nonzero=True)
directrix = Matrix([a, b, c])
focus = Matrix([fx, fy, 1])
# MatMul returns a scalar if none of its arguments is a matrix.
parabola = cast(
    "MatMul", MatMul(k, conic_from_focus_and_directrix(focus, directrix, 1))
)

println_indented(eq_chain(MatrixSymbol("P", 3, 3), parabola))

print("Its determinant:\n")

det = parabola.as_explicit().det().factor()
println_indented(eq_chain(Symbol("det"), det))

print("The determinant's factors are non-zero:\n")
print("  • n = a² + b² > 0 (the directrix is a finite line)")
print("  • d = a·fx + b·fy + c ≠ 0 (the focus is not on the directrix)\n")

n = Symbol("n", positive=True)
d = Symbol("d", nonzero=True)
det_sign = sign(det.subs({a * a + b * b: n, a * fx + b * fy + c: d}))
# SymPy doesn't simplify sign(k³) to sign(k) for real k.
det_sign_simple = det_sign.subs(sign(k**3), sign(k))
assert det_sign_simple == sign(k)

println_indented(eq_chain(sign(Symbol("det")), det_sign, det_sign_simple))

print(
    "Normalization multiplies a non-degenerate conic by its determinant sign.\n"
    "The upper left 2×2 block after normalization:\n"
)

# SymPy types Expr * Matrix as Expr (fix upstream).
block = cast("Matrix", det_sign_simple * parabola.as_explicit()[:2, :2])
block = block.subs(k * sign(k), Abs(k))
B = MatrixSymbol("B", 2, 2)
println_indented(eq_chain(B, block))

print("It's negative semidefinite:\n")

u, v = symbols("u v", real=True)
direction = Matrix([u, v])
block_form = quadratic_form(block, direction).factor()
println_indented(
    eq_chain(MatMul(direction.T, B, direction), Le(block_form, 0, evaluate=False))
)

assert block_form.is_nonpositive

print(HORIZONTAL_LINE + "\n")

print("Lemma: the normalization factor for conics with exactly one ideal point")
print("(parabolas, ideal point conics, parallel lines and double lines) is")
print("-sign(a + c), where M is the upper left 2×2 block of the conic matrix:\n")

block = Matrix([[a, b], [b, c]])
M = MatrixSymbol("M", 2, 2)
println_indented(eq_chain(M, block))

print("One ideal point: b² = a·c and M ≠ 0, hence a + c ≠ 0.")
print("Normalization makes M negative semidefinite: for parabolas by the lemma")
print("above, for the degenerate conics by definition. So does -sign(a + c):\n")

q = quadratic_form(block, direction)
# -(a+c)·q as sum of squares:
q_as_sum_squares = (
    -(a * q).expand().subs(a * c, b * b).factor()
    - (c * q).expand().subs(a * c, b * b).factor()
)

println_indented(
    eq_chain(
        MatMul(direction.T, -(a + c), M, direction, evaluate=False),
        Le(q_as_sum_squares, 0, evaluate=False),
    )
)
assert q_as_sum_squares.is_nonpositive
