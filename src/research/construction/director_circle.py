#!/usr/bin/env python

from sympy import simplify, sqrt, symbols

from lib.central_conic import conic_center, primary_radius, secondary_radius
from lib.circle import circle
from lib.matrix import conic_matrix
from research.sympy_utils import println_indented

conic = conic_matrix(*symbols("a,b,c,d,e,f"))

print("\nDirector circle radius:\n")

radius_square = primary_radius(conic) ** 2 + secondary_radius(conic) ** 2
radius = sqrt(radius_square.simplify().expand().factor(deep=True))
println_indented(radius.subs(conic.det(), symbols("det")))

print("\nDirector circle matrix:\n")

center = conic_center(conic)
director_circle = simplify(circle(center, r=radius))
println_indented(director_circle)

print("\nDirector circle matrix in adjugate form:\n")

a_adj, _, _, b_adj, c_adj, _, d_adj, e_adj, f_adj = conic.adjugate()
println_indented(
    director_circle.subs(
        zip(
            [a_adj, b_adj, c_adj, d_adj, e_adj, f_adj],
            symbols("a_adj b_adj c_adj d_adj e_adj f_adj"),
            strict=True,
        ),
    ),
)
