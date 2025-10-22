#!/usr/bin/env python

from sympy import pprint, symbols

from lib.intersection import conic_x_line
from lib.line import IDEAL_LINE
from lib.matrix import conic_matrix

print("\nIdeal points on a conic:\n")

pprint(conic_x_line(conic_matrix(*symbols("a,b,c,d,e,f")), IDEAL_LINE))

print("\nPotential formulae based on the coefficients' signs:\n")

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
            conic = conic_matrix(a, b, c, d, e, f)
            ideal_points = conic_x_line(conic, IDEAL_LINE)
            if str(ideal_points) not in solutions:
                solutions.add(str(ideal_points))
                print()
                pprint(ideal_points)
