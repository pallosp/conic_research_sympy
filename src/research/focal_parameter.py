#!/usr/bin/env python

from sympy import sqrt, symbols

from lib.distance import PointLineDistance
from lib.matrix import ConicMatrix
from lib.parabola import ParabolaDirectrix, ParabolaFocus
from research.util import println_indented

print("\nParabola focal parameter (focus-directrix distance):\n")

a, b, c, d, e, f = symbols("a b c d e f", real=True)
parabola = ConicMatrix(a, b, c, d, e, f)
det = parabola.det().subs(b * b, a * c)

focal_parameter = (
    PointLineDistance(ParabolaFocus(parabola), ParabolaDirectrix(parabola))
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
