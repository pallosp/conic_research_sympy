#!/usr/bin/env python

from typing import cast

from sympy import Abs, Expr, Symbol, expand, sign, sqrt, symbols

from lib.central_conic import linear_eccentricity, primary_radius
from lib.conic import projective_conic_center
from lib.conic_direction import ConicNormFactor
from lib.distance import point_line_distance
from lib.matrix import conic_matrix
from lib.parabola import parabola_directrix, parabola_focus
from lib.type_aliases import Conic
from research.sympy_utils import println_indented

HORIZONTAL_LINE = "-" * 80

print("\nParabola focal parameter (focus-directrix distance):\n")

a, b, c, d, e, f = symbols("a b c d e f", real=True)
det = Symbol("det", nonzero=True)
parabola = conic_matrix(a, b, c, d, e, f)
parabola_det = parabola.det().subs(b * b, a * c)

parabola_fp = (
    point_line_distance(parabola_focus(parabola), parabola_directrix(parabola))
    .subs(b * b, a * c)
    .factor()
    .subs(b * b, a * c)
    .factor()
)

# Take the absolute value
parabola_fp = sqrt(parabola_fp**2)

println_indented(parabola_fp.subs(parabola_det, det))

print("Alternative form:\n")

cx, cy, _ = projective_conic_center(parabola).flat()
parabola_fp = parabola_fp.subs(parabola_det, -(cx**2 + cy**2) / (a + c)).simplify()

println_indented(parabola_fp)

print(HORIZONTAL_LINE + "\n")


def _simplify_focal_param(conic: Conic, focal_param: Expr, abs_disc: Expr) -> Expr:
    """Simplifies a central conic's focal parameter formula.

    `conic` must be built from the module-level `a`, ..., `f` symbols.
    `abs_disc` is the value of `|ac - b²|`.
    """
    disc_factorized = (
        (a + c + sign(det) * sqrt((a - c) ** 2 + 4 * b**2))
        * (a + c - sign(det) * sqrt((a - c) ** 2 + 4 * b**2))
        / 4
    )
    assert disc_factorized.expand() == a * c - b * b

    focal_param = (
        focal_param.subs(conic.det(), det)
        .subs(ConicNormFactor(conic), sign(det))
        .subs(Abs(a * c - b * b), abs_disc)
        .subs(a * c - b * b, disc_factorized)
        .subs(Abs(det), det * sign(det))
        .simplify()
    )
    return cast("Expr", focal_param)


def hyperbola_focal_param(conic: Conic) -> Expr:
    """Computes a hyperbola's focal parameter.

    `conic` must be built from the module-level `a`, ..., `f` symbols.

    hyperbola point order:
      center | directrix | vertex | focus
      C ... D ... V ... F

    linecc = F - C
    r = V - C
    ecc = (F - V) / (V - D)
    ecc = linecc / r

    F - V = linecc - r
    V - D = (linecc - r) / ecc

    focal parameter
      = F - D
      = (F - V) + (V - D)
      = (linecc - r) * (1 + 1 / ecc)
      = (linecc - r) * (1 + r / linecc)
      = (linecc**2 - r**2) / linecc
    """
    linecc = linear_eccentricity(conic)
    r = primary_radius(conic)
    focal_param = (linecc**2 - r**2) / linecc
    return _simplify_focal_param(conic, focal_param, b * b - a * c)


def ellipse_focal_param(conic: Conic) -> Expr:
    """Computes an ellipse's focal parameter.

    `conic` must be built from the module-level `a`, ..., `f` symbols.

    ellipse point order:
      center | focus | vertex | directrix
      C ... F ... V ... D

    linecc = F - C
    r = V - C
    ecc = (V - F) / (D - V)
    ecc = linecc / r

    V - F = r - linecc
    D - V = (r - linecc) / ecc

    focal parameter
      = D - F
      = (D - V) + (V - F)
      = (r - linecc) * (1 + 1 / ecc)
      = (r - linecc) * (1 + r / linecc)
      = (r² - linecc²) / linecc
    """
    linecc = linear_eccentricity(conic)
    r = primary_radius(conic)
    focal_param = (r**2 - linecc**2) / linecc
    return _simplify_focal_param(conic, focal_param, a * c - b * b)


generic_conic = conic_matrix(a, b, c, d, e, f)

print("Hyperbola focal parameter:\n")

hyperbola_formula = hyperbola_focal_param(generic_conic)
println_indented(hyperbola_formula)

print("Ellipse focal parameter:\n")

ellipse_formula = ellipse_focal_param(generic_conic)
println_indented(ellipse_formula)

assert ellipse_formula == hyperbola_formula
print("The two formulas are identical.\n")

print("They even work for parabolas:\n")

parabola_formula = (
    hyperbola_formula.subs(b * b, a * c)
    .subs(det, -(cx**2 + cy**2) / (a + c))
    .factor(deep=True)
    # Identity that steers simplify() towards cancelling the |a+c| factors
    .subs(Abs(a + c), (a + c) ** 2 / Abs(a + c))
    .simplify()
    .factor()
    # -(a+c)*det = cx²+cy² ≠ 0
    #  => a+c ≠ 0
    #  => |√(a+c) * √(1/(a+c))| = 1
    .subs(Abs(sqrt(a + c) * sqrt(1 / (a + c))), 1)
    .subs(expand(cx**2), cx**2)
    # Same sign convention as in the alternative form above
    .subs(expand(cy**2), (-cy) ** 2)
)

println_indented(parabola_formula)

assert parabola_formula == parabola_fp
