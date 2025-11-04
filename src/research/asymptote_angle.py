#!/usr/bin/env python

from sympy import Eq, Limit, cos, oo, solve, symbols

from lib.sympy_utils import eq_chain
from research.util import println_indented

print("\nAngle between the focal axis and an asymptote of a hyperbola\n")

print(
    "Take the triangle among the two foci (f₁, f₂) of a hyperbola "
    "and a point (p) on it:\n",
)

print("  distance(f₁, f₂) = 2 * linear_eccentricity")
print("  distance(f₁, p) = x")
print("  distance(f₂, p) = 2 * primary_radius + x")

print("\nLaw of cosines to compute the angle at f₁:\n")

lin_ecc, radius, x, cos_f1 = symbols("l r x cos(f₁)")
law_of_cosines = Eq(
    (2 * radius + x) ** 2,
    x**2 + (2 * lin_ecc) ** 2 - 2 * x * (2 * lin_ecc) * cos_f1,
)
println_indented(law_of_cosines)

solution = solve(law_of_cosines, cos_f1)[0]
println_indented(Eq(cos_f1, solution))

print("The angle between the focal axis and an asymptote:\n")

cos_f1_limit = Limit(cos_f1, x, oo)
aa_angle = symbols("axis_asymptote_angle")
println_indented(
    eq_chain(
        cos(aa_angle),
        -cos_f1_limit,
        -cos_f1_limit.subs(cos_f1, solution),
        -cos_f1_limit.subs(cos_f1, solution).doit(),
    ),
)
