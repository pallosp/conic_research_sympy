#!/usr/bin/env python

from sympy import pprint, symbols

from lib.intersection import ConicXLine
from lib.line import IDEAL_LINE
from lib.matrix import ConicMatrix

print("\nIdeal points on a conic:")

solutions = set()
for a_assumption in [{"zero": True}, {"nonzero": True}]:
    for b_assumption in [{"negative": True}, {"zero": True}, {"positive": True}]:
        for c_assumption in [{"zero": True}, {"nonzero": True}]:
            a = symbols("a", **a_assumption)
            b = symbols("b", **b_assumption)
            c = symbols("c", **c_assumption)
            if a.equals(0) and b.equals(0) and c.equals(0):
                continue
            d, e, f = symbols("d,e,f")
            conic = ConicMatrix(a, b, c, d, e, f)
            ideal_points = ConicXLine(conic, IDEAL_LINE)
            if str(ideal_points) not in solutions:
                solutions.add(str(ideal_points))
                print()
                pprint(ideal_points)
