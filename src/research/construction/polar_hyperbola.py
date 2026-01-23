#!/usr/bin/env python

from sympy import I, cos, factor, simplify, sin, symbols

from lib.transform import homography_from_samples
from research.sympy_utils import println_indented

print()
print(
    "Hyperbola in polar form, expressed with the center, "
    "the focal axis angle and the radii:"
)
print()

cx, cy, r1, r2, alpha = symbols("cx cy r1 r2 alpha")

vertex_1 = (cx + r1 * cos(alpha), cy + r1 * sin(alpha))
vertex_2 = (cx - r1 * cos(alpha), cy - r1 * sin(alpha))

# Formula: conic_properties/ideal_points.py
ideal_point_1 = (
    r1 * cos(alpha) - I * r2 * sin(alpha),
    r1 * sin(alpha) + I * r2 * cos(alpha),
    0,
)
ideal_point_2 = (
    r1 * cos(alpha) + I * r2 * sin(alpha),
    r1 * sin(alpha) - I * r2 * cos(alpha),
    0,
)

source = ((1, 0), (0, 1), (-1, 0), (0, -1))
target = (vertex_1, ideal_point_1, vertex_2, ideal_point_2)

hyperbola = homography_from_samples(source, target)
hyperbola = hyperbola.applyfunc(factor).applyfunc(simplify)
hyperbola /= hyperbola[6]

println_indented(hyperbola)
