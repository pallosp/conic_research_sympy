#!/usr/bin/env python

from sympy import pprint, symbols

from lib.distance import PointLineDistance
from lib.matrix import ConicMatrix
from lib.parabola import ParabolaDirectrix, ParabolaFocus

print("\nParabola focal parameter (focus-directrix distance):\n")

a, b, c, d, e, f = symbols("a b c d e f")
parabola = ConicMatrix(a, b, c, d, e, f)
det = parabola.det().subs(b * b, a * c)

focal_parameter = (
    PointLineDistance(ParabolaFocus(parabola), ParabolaDirectrix(parabola))
    .subs(b * b, a * c)
    .factor()
    .subs(b * b, a * c)
    .factor()
    .subs(det, symbols("det"))
)

pprint(focal_parameter)
