#!/usr/bin/env python

from sympy import expand, factor, Matrix, pprint, simplify, sqrt, symbols

from lib.circle import Circle
from lib.central_conic import ConicCenter, SemiMajorAxis, SemiMinorAxis


a, b, c, d, e, f = symbols("a,b,c,d,e,f")
conic = Matrix([[a, b, d], [b, c, e], [d, e, f]])

center_x, center_y = ConicCenter(conic)
radius_square = SemiMajorAxis(conic) ** 2 + SemiMinorAxis(conic) ** 2
radius_square = factor(expand(simplify(radius_square)))
radius = sqrt(radius_square)

director_circle = simplify(Circle(center_x, center_y, radius))

print("\nDirector circle radius:")
pprint(radius)
print("\nDirector circle matrix:")
pprint(director_circle)
