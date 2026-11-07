#!/usr/bin/env python

"""Derives the latus rectum formula from the eccentricity and the focal parameter."""

from sympy import Abs, Eq, Symbol, expand, factor, radsimp, symbols

from lib.central_conic import primary_radius, secondary_radius
from lib.circle import circle_radius
from lib.conic import eccentricity, focal_parameter
from lib.conic_direction import ConicNormFactor
from lib.matrix import conic_matrix
from research.sympy_utils import eq_chain, println_indented

a, b, c, d, e, f = symbols("a:f", real=True)
conic = conic_matrix(a, b, c, d, e, f)
lr_symbol = Symbol("latus rectum")
det_symbol = Symbol("det", real=True)

print("\nLatus rectum of non-circular conics:\n")
println_indented(Eq(lr_symbol, 2 * Symbol("eccentricity") * Symbol("focal param")))

fp = focal_parameter(conic)
ecc = eccentricity(conic)
lr = 2 * ecc * fp

lr = lr.subs(conic.det(), det_symbol)
lr = lr.subs(Abs(det_symbol), det_symbol * ConicNormFactor(conic))
lr_squared = (lr**2).factor()

# Rewrite the fraction to a form that sympy can simplify.
norm_sign = ConicNormFactor(conic)
lr_squared = lr_squared.replace(
    lambda x: x.is_Add and x.has(norm_sign),
    lambda x: norm_sign * expand(norm_sign * x),
)

println_indented(Eq(Symbol("latus rectum") ** 2, lr_squared))

print(
    "The simplified formula also works for circles, "
    "where eccentricity = 0 and focal param = ∞.\n"
)

lr_squared_circle = lr_squared.subs(b, 0).subs(c, a)

println_indented(
    eq_chain(Symbol("latus rectum") ** 2, lr_squared_circle, Symbol("diameter") ** 2)
)

diameter = 2 * circle_radius(conic)
assert lr_squared_circle == diameter.subs(conic.det(), det_symbol) ** 2

print("Alternative, radius-based formula for central conics:\n")
println_indented(Eq(lr_symbol, 2 * Symbol("r_secondary") ** 2 / Symbol("r_primary")))

lr_alt = 2 * (secondary_radius(conic) ** 2 / primary_radius(conic))

lr_alt_squared = lr_alt**2
lr_alt_squared = factor(1 / radsimp(1 / lr_alt_squared))
lr_alt_squared = lr_alt_squared.subs(conic.det(), det_symbol)

println_indented(Eq(lr_symbol**2, lr_alt_squared))

assert lr_alt_squared == lr_squared
