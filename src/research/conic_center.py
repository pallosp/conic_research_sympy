#!/usr/bin/env python3

from sympy import Eq, Matrix, Poly, pprint, solve, symbols

from lib.matrix import ConicMatrix, QuadraticForm

x, y, cx, cy = symbols("x y cx cy", real=True)

conic_matrix = ConicMatrix(*symbols("a b c d e f"))
conic = Poly(QuadraticForm(conic_matrix, Matrix([x, y, 1])), [x, y])

# Reflect conic to point (cx, cy)
reflected = Poly(
    QuadraticForm(conic_matrix, Matrix([2 * cx - x, 2 * cy - y, 1])), [x, y],
)

# Quadratic part must be invariant under reflection
for mon in (x**2, x * y, y**2):
    lhs, rhs = conic.coeff_monomial(mon), reflected.coeff_monomial(mon)
    assert lhs == rhs, f"Quadratic mismatch for {mon}: {lhs} != {rhs}"

# Solve for (cx, cy) by matching linear terms
solution = solve(
    [
        Eq(conic.coeff_monomial(x), reflected.coeff_monomial(x)),
        Eq(conic.coeff_monomial(y), reflected.coeff_monomial(y)),
    ],
    (cx, cy),
)

# Verify constant term consistency
orig_const = conic.coeff_monomial(1)
refl_const = reflected.coeff_monomial(1).subs(solution).simplify()
assert orig_const == refl_const

# Result
print("\nConic center point:\n")
pprint(solution)
