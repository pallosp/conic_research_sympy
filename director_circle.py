#!/usr/bin/env python

from sympy import Matrix, pprint, simplify, sqrt, symbols

from lib.circle import Circle
from lib.central_conic import ConicCenter, SemiMajorAxis, SemiMinorAxis


a, b, c, d, e, f = symbols("a,b,c,d,e,f")
conic = Matrix([[a, b, d], [b, c, e], [d, e, f]])
center_x, center_y = ConicCenter(conic)
radius = sqrt(SemiMajorAxis(conic) ** 2 + SemiMinorAxis(conic) ** 2)
director_circle = Circle(center_x, center_y, radius)
pprint(simplify(director_circle))
