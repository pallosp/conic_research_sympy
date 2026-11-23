#!/usr/bin/env python

from sympy import Abs, Eq, Symbol, sqrt, symbols

from lib.conic import eccentricity, focal_parameter
from lib.conic_direction import ConicNormFactor
from lib.matrix import conic_matrix
from lib.parabola import parabola_focal_length
from research.sympy_utils import println_indented

a, b, c, d, e, f = symbols("a:f", real=True)
conic = conic_matrix(a, b, c, d, e, f)

print("\nConic focal length (focus-vertex distance):\n")

fp = focal_parameter(conic)
ecc = eccentricity(conic)
fl = fp * ecc / (ecc + 1)

det = Symbol("det", real=True)
norm_sign = Symbol("norm_sign", nonzero=True)
s = Symbol("s", nonnegative=True)
q = Symbol("q", real=True)

fl = (
    fl.subs(ConicNormFactor(conic), norm_sign)
    .subs(((a - c) ** 2 + 4 * b**2).expand(), s)
    .subs((a - c) ** 2 + 4 * b**2, s)
    .subs(conic.det(), det)
    .subs(sqrt(s) - norm_sign * (a + c), q)
    .subs(sqrt(1 / q), sqrt(q) / Abs(q))
    .expand()
    .subs(q**2 / Abs(q), Abs(q))
)

println_indented(Eq(s, (a - c) ** 2 + 4 * b**2))
println_indented(Eq(q, sqrt(s) - ConicNormFactor(conic) * (a + c)))
println_indented(Eq(Symbol("focal_length"), fl))

print("Simplified formula for parabolas:\n")

# b² = ac  =>  s = (a + c)²
# ConicNormFactor = -sign(a + c)  =>  q = 2|a + c|
abs_apc = Symbol("|a + c|", positive=True)
fl = fl.subs(s, abs_apc**2).subs(q, 2 * abs_apc)

println_indented(Eq(Symbol("focal_length"), fl))

# Cross-check with parabola_focal_length() using sign(det) = -sign(a + c)
apc_sign = Symbol("apc_sign", nonzero=True)
parabola_fl = parabola_focal_length(conic).subs(conic.det(), det)
assert (
    (fl**2 - parabola_fl**2)
    .subs(Abs(det), -det * apc_sign)
    .subs(abs_apc, apc_sign * (a + c))
    .subs(apc_sign**2, 1)
    .simplify()
) == 0
