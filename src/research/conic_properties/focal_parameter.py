#!/usr/bin/env python

from sympy import sqrt, symbols

from lib.distance import point_line_distance
from lib.matrix import conic_matrix
from lib.parabola import parabola_directrix, parabola_focus
from research.sympy_utils import println_indented

print("\nParabola focal parameter (focus-directrix distance):\n")

a, b, c, d, e, f = symbols("a b c d e f", real=True)
parabola = conic_matrix(a, b, c, d, e, f)
det = parabola.det().subs(b * b, a * c)

focal_parameter = (
    point_line_distance(parabola_focus(parabola), parabola_directrix(parabola))
    .subs(b * b, a * c)
    .factor()
    .subs(b * b, a * c)
    .factor()
)

# Take the absolute value
focal_parameter = sqrt(focal_parameter**2)

println_indented(focal_parameter.subs(det, symbols("det")))

print("Alternative form:\n")

cx, cy, _ = parabola.row(0).cross(parabola.row(1))
focal_parameter = focal_parameter.subs(det, -(cx**2 + cy**2) / (a + c)).simplify()

println_indented(focal_parameter)
