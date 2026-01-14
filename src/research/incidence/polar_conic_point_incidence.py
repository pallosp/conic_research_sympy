#!/usr/bin/env python3

from sympy import Matrix, symbols

from lib.matrix import quadratic_form
from lib.polar_conic import conic_from_polar_matrix
from research.sympy_utils import println_indented

# First approach:
#
#   P [cos α, sin α, 1]ᵀ = [x, y, z]ᵀ
#   λ [cos α, sin α, 1]ᵀ = P⁻¹ [x, y, z]ᵀ
#
# When P⁻¹ [x, y, z]ᵀ = [ls, lc, l1]ᵀ, there must be a solution for α, which exists
# iff ls² + lc² = l1².

polar_conic = Matrix(3, 3, symbols("a b c d e f g h i"))
point = Matrix(symbols("x y z"))
l_times_cos_a, l_times_sin_a, l_times_1 = polar_conic.adjugate() * point
incidence_poly_1 = l_times_1**2 - (l_times_cos_a**2 + l_times_sin_a**2)

# Second approach:
#
# Convert the polar conic to cartesian form, then use the quadratic form to
# check point incidence.

conic = conic_from_polar_matrix(polar_conic)
incidence_poly_2 = quadratic_form(conic, point)

# The two approaches are equivalent.

assert incidence_poly_1.expand() == incidence_poly_2.expand()

print("\nA polar conic contains point (x,y,z) iff the following expression is zero:\n")

println_indented(incidence_poly_1)

# The incidence polynomial can't be factorized.
