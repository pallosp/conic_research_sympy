#!/usr/bin/env python

from sympy import Abs, Eq, I, Limit, cos, oo, sign, sin, solve, sqrt, symbols, tan

from lib.central_conic import primary_radius, secondary_radius
from lib.conic_direction import ConicNormFactor
from lib.matrix import conic_matrix
from research.sympy_utils import eq_chain, println_indented

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

print("Angle between the focal axis and an asymptote:\n")

# cos

cos_f1_limit = Limit(cos_f1, x, oo)
aa_angle = symbols("axis_asymptote_angle")
println_indented(
    eq_chain(
        cos(aa_angle),
        -cos_f1_limit,
        -cos_f1_limit.subs(cos_f1, solution),
        -cos_f1_limit.subs(cos_f1, solution).doit(),
        1 / symbols("eccentricity"),
    ),
)

# sin

sec_radius = symbols("r_2")
println_indented(
    eq_chain(
        sin(aa_angle),
        sqrt(lin_ecc**2 - radius**2) / lin_ecc,
        sqrt(-(sec_radius**2)) / lin_ecc,
        sec_radius * I / lin_ecc,
    )
)

# tan

a, b, c, d, e, f = symbols("a b c d e f", real=True)
conic = conic_matrix(a, b, c, d, e, f)
det = symbols("det", nonzero=True)
tan_aa = I * secondary_radius(conic) / primary_radius(conic)
tan_aa = sqrt(tan_aa**2).simplify().subs(ConicNormFactor(conic), sign(det))
num, denom = (tan_aa**2).as_numer_denom()
println_indented(
    eq_chain(
        tan(aa_angle),
        I * sec_radius / radius,
        tan_aa,
        Abs(num) / sqrt(num * denom).expand().factor(),
    )
)
