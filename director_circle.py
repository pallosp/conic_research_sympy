#!/usr/bin/env python

from sympy import expand, factor, pprint, simplify, sqrt, symbols

from lib.circle import Circle
from lib.central_conic import ConicCenter, SemiMajorAxis, SemiMinorAxis
from lib.matrix import ConicMatrix


conic = ConicMatrix(*symbols("a,b,c,d,e,f"))

center = ConicCenter(conic)
radius_square = SemiMajorAxis(conic) ** 2 + SemiMinorAxis(conic) ** 2
radius_square = factor(expand(simplify(radius_square)))
radius = sqrt(radius_square)

director_circle = simplify(Circle(center, radius))

print("\nDirector circle radius:")
pprint(radius)
print("\nDirector circle matrix:")
pprint(director_circle)

print("\nDirector circle matrix in adjugate form:")
a_adj, b_adj, d_adj, _, c_adj, e_adj, _, _, f_adj = conic.adjugate()
pprint(
    director_circle.subs(
        zip(
            [a_adj, b_adj, c_adj, d_adj, e_adj, f_adj],
            symbols("a_adj,b_adj,c_adj,d_adj,e_adj,f_adj"),
        ),
    )
)
