#!/usr/bin/env python

from sympy import Abs, Eq, pprint, symbols

from lib.central_conic import linear_eccentricity, primary_radius, secondary_radius
from lib.conic_direction import ConicNormFactor
from lib.matrix import conic_matrix

conic = conic_matrix(*symbols("a b c d e f", real=True))
det = symbols("det", real=True)

le_square = (linear_eccentricity(conic) ** 2).subs(conic.det(), det).factor()

pprint(Eq(symbols("linear_eccentricity") ** 2, le_square))
print()

r1 = primary_radius(conic)
r2 = secondary_radius(conic)

radius_square_diff = (
    (r1**2 - r2**2)
    .simplify()
    .factor()
    .subs(conic.det(), det)
    .subs(ConicNormFactor(conic) * det, Abs(det))
)

pprint(Eq(symbols("r1") ** 2 - symbols("r2") ** 2, radius_square_diff))
print()
