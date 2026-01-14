#!/usr/bin/env python

from sympy import Matrix, pprint, simplify, symbols

from lib.line import line_through_point
from lib.point import ideal_point_on_line

line1 = Matrix(symbols("a1,b1,c1"))
line2 = Matrix(symbols("a2,b2,c2"))

vertex = line1.cross(line2)

# Negate one of the lines in order to choose the bisector whose points
# substituted into the lines' equations are positive.
dir1 = ideal_point_on_line(line1)
dir2 = ideal_point_on_line(-line2)
dir_bisector = dir1.normalized() + dir2.normalized()

angle_bisector = line_through_point(vertex, direction=dir_bisector)
angle_bisector = simplify(angle_bisector)

print("\nAngle bisector of two lines:\n")
pprint(angle_bisector)
