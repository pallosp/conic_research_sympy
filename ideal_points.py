#!/usr/bin/env python

from sympy import expand, factor, pprint, simplify, symbols, Matrix
from sympy.solvers import solve

from lib.intersection import ConicXLine
from lib.line import IDEAL_LINE
from lib.matrix import ConicMatrix, QuadraticForm

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
            ideal_points_str = "{0},{1}".format(*ideal_points)
            if ideal_points_str not in solutions:
                solutions.add(ideal_points_str)
                print()
                pprint(ideal_points)
