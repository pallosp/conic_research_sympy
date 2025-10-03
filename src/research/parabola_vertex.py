#!/usr/bin/env python

from sympy import Matrix, pprint, simplify, sqrt, symbols

from lib.matrix import ConicMatrix

print("\nParabola vertex:\n")

a, b, c, d, e, f = symbols("a b c d e f", real=True)
parabola = ConicMatrix(a, b, c, d, e, f)

focus = Matrix(symbols("focus.x focus.y"))

d_a, e_a, _ = parabola.row(0).cross(parabola.row(1))
focal_parameter = sqrt(d_a**2 + e_a**2) / (a + c) ** 2

axis_direction = Matrix([d_a, e_a])
normalized_axis_direction = axis_direction / axis_direction.norm()

vertex = focus - (normalized_axis_direction * focal_parameter / 2).applyfunc(simplify)

pprint(vertex)
