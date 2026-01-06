#!/usr/bin/env python

from sympy import expand, factor, pprint, simplify, symbols

from lib.circle import UNIT_CIRCLE
from lib.transform import rotate, scale_xy, transform_conic, translate

center = symbols("cx,cy")
r1, r2, angle = symbols("r1 r2 theta")
scaling = scale_xy(r1, r2)
rotation = rotate(angle)
translation = translate(center)
transformation = translation * rotation * scaling
ellipse = transform_conic(UNIT_CIRCLE, transformation)

a, b, c = symbols("a,b,c")
for i in (2, 5, 6, 7, 8):
    ellipse[i] = factor(ellipse[i], *center)
    ellipse[i] = ellipse[i].subs(
        [
            (ellipse[0], "a"),
            (expand(ellipse[1]), "b"),
            (ellipse[4], "c"),
        ],
    )
    ellipse[i] = factor(simplify(ellipse[i]), *center)
    ellipse[i] = ellipse[i].subs(expand(simplify(ellipse[1]) * 2), "b*2")

print("\nEllipse(cx, cy, r1, r2, θ):\n")
pprint(ellipse)
