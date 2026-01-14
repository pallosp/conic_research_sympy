#!/usr/bin/env python

from sympy import pprint, simplify, sqrt, symbols

from lib.central_conic import conic_center, primary_radius, secondary_radius
from lib.circle import circle
from lib.matrix import conic_matrix

conic = conic_matrix(*symbols("a,b,c,d,e,f"))

center = conic_center(conic)
radius_square = primary_radius(conic) ** 2 + secondary_radius(conic) ** 2
radius = sqrt(radius_square.simplify().factor())

print("\nDirector circle radius:")
det = symbols("det")
pprint(radius.subs(conic.det(), det))

print("\nDirector circle matrix:")
director_circle = simplify(circle(center, radius))
pprint(director_circle)

print("\nDirector circle matrix in adjugate form:")
a_adj, b_adj, d_adj, _, c_adj, e_adj, _, _, f_adj = conic.adjugate()
pprint(
    director_circle.subs(
        zip(
            [a_adj, b_adj, c_adj, d_adj, e_adj, f_adj],
            symbols("a_adj,b_adj,c_adj,d_adj,e_adj,f_adj"),
            strict=False,
        ),
    ),
)
