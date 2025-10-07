#!/usr/bin/env python

from sympy import Eq, simplify, solve, sqrt, symbols

from lib.central_conic import ConicCenter
from lib.matrix import ConicMatrix
from lib.transform import Rotate, TransformConic, Translate
from research.util import println_indented

# Bring the conic to the standard form by translating its center to the origin
# and rotating it to align with the axes.

a, b, c, d, e, f = symbols("a b c d e f")
conic = ConicMatrix(a, b, c, d, e, f)
det = conic.det()

# Translate the center to the origin
cx, cy = ConicCenter(conic)
conic = simplify(TransformConic(conic, Translate(-cx, -cy)))

# Rotate with theta
theta = symbols("theta")
conic = simplify(TransformConic(conic, Rotate(theta)))

# Choose an angle for which the conic will be axis-aligned, i.e. the xy term
# vanishes.
angle = solve(conic[1], theta)[1]

# Extract √((a-c)² + 4b²) to a variable to help sympy simplify the result
subexpr_val = sqrt((a - c) ** 2 + 4 * b**2).expand()
subexpr_var = symbols("s")
angle = angle.subs(subexpr_val, subexpr_var)

conic = conic.subs(det, symbols("det"))
conic = conic.subs(theta, angle).simplify()
conic = conic.subs(4 * b**2, 4 * b**2 + subexpr_var**2 - subexpr_val**2).simplify()
conic = conic.subs(subexpr_var, subexpr_val)

# After the rotation the conic should be diagonal.
assert conic.is_diagonal()

# Standard central conic formula: x²/r1² + y²/r2² = 1
print("\nSemi-axis lengths:\n")
a, c, f = conic.diagonal()
println_indented(Eq(symbols("r1") ** 2, -f / a))
println_indented(Eq(symbols("r2") ** 2, -f / c))

print("In terms of the 2x2 submatrix eigenvalues:\n")

eigenvals = list(conic[:2, :2].eigenvals())
eigenval_subs = {
    eigenvals[0]: symbols("lambda_1"),
    eigenvals[1]: symbols("lambda_2"),
}

println_indented(Eq(symbols("r1") ** 2, -f / a.subs(eigenval_subs)))
println_indented(Eq(symbols("r2") ** 2, -f / c.subs(eigenval_subs)))
