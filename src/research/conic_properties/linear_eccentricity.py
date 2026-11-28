#!/usr/bin/env python

from sympy import Abs, Symbol, pprint, symbols

from lib.central_conic import linear_eccentricity, primary_radius, secondary_radius
from lib.conic import eccentricity
from lib.conic_direction import ConicNormFactor
from lib.matrix import conic_matrix
from research.equations import eq_chain
from research.pretty_print import HORIZONTAL_LINE

conic = conic_matrix(*symbols("a:f", real=True))
det_symbol = Symbol("det", real=True)
linecc_symbol = Symbol("linear_eccentricity", real=True)
ecc_symbol = Symbol("eccentricity", real=True)
r1_symbol = Symbol("r1")
r2_symbol = Symbol("r2")

# Proof: linear_eccentricity² = r₁² - r₂²

le_square = (linear_eccentricity(conic) ** 2).subs(conic.det(), det_symbol).factor()

r1 = primary_radius(conic)
r2 = secondary_radius(conic)

radius_square_diff = (
    (r1**2 - r2**2)
    .simplify()
    .factor()
    .subs(conic.det(), det_symbol)
    .subs(ConicNormFactor(conic) * det_symbol, Abs(det_symbol))
)

assert le_square == radius_square_diff

print()
pprint(
    eq_chain(
        linecc_symbol**2,
        le_square,
        r1_symbol**2 - r2_symbol**2,
    )
)

print()
print(HORIZONTAL_LINE)
print()

# Proof: linear_eccentricity² = eccentricity² · r₁²

ecc_square_x_radius_square = (
    (eccentricity(conic) ** 2 * r1**2)
    .simplify()
    .factor()
    .subs(conic.det(), det_symbol)
    .subs(ConicNormFactor(conic) * det_symbol, Abs(det_symbol))
)

assert le_square == ecc_square_x_radius_square

pprint(
    eq_chain(
        linecc_symbol**2,
        ecc_square_x_radius_square,
        ecc_symbol**2 * r1_symbol**2,
    )
)
print()
