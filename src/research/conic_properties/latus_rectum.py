#!/usr/bin/env python

"""Derives the latus rectum formula in two ways: from the eccentricity and the
focal parameter, and from the radii of central conics.
"""

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

# Claim: semi-latus rectum = eccentricity · focal parameter
#
# Proof:
#  - F: a focus
#  - L: an endpoint of the latus rectum through F
#  - F', L': the projections of F and L onto the directrix corresponding to F
#
# FL is parallel to the directrix, so FLL'F' is a rectangle, and
#
#   |LL'| = |FF'| = focal parameter
#
# L lies on the conic, so |FL| / |LL'| = eccentricity. Therefore
#
#   semi-latus rectum = |FL| = eccentricity · |LL'| = eccentricity · focal parameter
#
# The formula breaks down for circles, where eccentricity = 0 and
# focal parameter = ∞.

println_indented(Eq(lr_symbol, 2 * Symbol("eccentricity") * Symbol("focal param")))

fp = focal_parameter(conic)
ecc = eccentricity(conic)
lr = 2 * ecc * fp

lr = lr.subs(conic.det(), det_symbol)
lr = lr.subs(Abs(det_symbol), det_symbol * ConicNormFactor(conic))
lr_squared = (lr**2).factor()

# Rewrite each sum X containing norm_sign as norm_sign·(norm_sign·X), which
# holds because norm_sign² = 1. SymPy then reduces norm_sign³ to norm_sign and
# cancels it.
norm_sign = ConicNormFactor(conic)
lr_squared = lr_squared.replace(
    lambda x: x.is_Add and x.has(norm_sign),
    lambda x: norm_sign * expand(norm_sign * x),
)

println_indented(Eq(lr_symbol**2, lr_squared))

print(
    "The simplified formula also works for circles, "
    "where eccentricity = 0 and focal parameter = ∞.\n"
)

lr_squared_circle = lr_squared.subs(b, 0).subs(c, a)

println_indented(eq_chain(lr_symbol**2, lr_squared_circle, Symbol("diameter") ** 2))

diameter = 2 * circle_radius(conic)
assert lr_squared_circle == diameter.subs(conic.det(), det_symbol) ** 2

print("Alternative, radius-based formula for central conics:\n")

# Claim: semi-latus rectum = r_secondary² / r_primary
#
# Proof for ellipses:
#  - O: the center
#  - F: a focus
#  - C: a co-vertex
#  - C': where the ray OC meets the auxiliary circle
#  - L: the endpoint of the latus rectum through F, on the same side as C
#  - L': the fourth vertex of the rectangle OFL'C
#
#                    ......C'.....
#               .....      :      .....
#             ..           :           ..
#           ..             :             ..
#         ..       ********C********-------L'
#        .    *****        |        *****  | .
#       .  ***             |             **L  .
#       .**                |               |**.
#      .*                  |               |  *.
#      *                   |               |   *
#      *-------------------O---------------F---*
#
#    * ellipse    . auxiliary circle
#
#  1. |OC| = r_secondary and |OC'| = r_primary.
#  2. C is equally far from the two foci, and the distances add up to
#     2·r_primary, so |FC| = r_primary. The diagonals of a rectangle are equal,
#     so |OL'| = |FC| = r_primary, i.e. L' lies on the auxiliary circle.
#  3. Stretching the plane perpendicular to the major axis by the factor
#     r_primary / r_secondary maps the ellipse to its auxiliary circle. It fixes
#     O and F, maps C to C', and maps L to L'. It scales every segment
#     perpendicular to the major axis by the same factor, so
#     |FL| / |FL'| = |OC| / |OC'|.
#  4. |FL'| = |OC|, as opposite sides of the rectangle. Therefore
#
#       semi-latus rectum = |FL| = |OC|² / |OC'| = r_secondary² / r_primary
#
# For hyperbolas, r_secondary is imaginary. The algebraic check below covers
# them too.

println_indented(Eq(lr_symbol, 2 * Symbol("r_secondary") ** 2 / Symbol("r_primary")))

lr_alt = 2 * (secondary_radius(conic) ** 2 / primary_radius(conic))

# Rationalize the numerator: radsimp multiplies the reciprocal's denominator by
# its conjugate, turning it into a multiple of a·c - b², which then cancels.
lr_alt_squared = lr_alt**2
lr_alt_squared = factor(1 / radsimp(1 / lr_alt_squared))
lr_alt_squared = lr_alt_squared.subs(conic.det(), det_symbol)

println_indented(Eq(lr_symbol**2, lr_alt_squared))

assert lr_alt_squared == lr_squared
