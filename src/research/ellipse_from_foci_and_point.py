#!/usr/bin/env python

from sympy import cancel, expand, gcd, sqrt, symbols

from lib.central_conic import conic_from_foci_and_radius
from lib.distance import point_point_distance
from lib.ellipse import ellipse
from lib.point import centroid
from research.sympy_utils import println_indented

print("\nEllipse from its focus points and an incident point:\n")

cx, cy = symbols("cx cy")  # center
dx, dy = symbols("dx dy")  # center -> focus vector

focus1 = (cx + dx, cy + dy)
focus2 = (cx - dx, cy - dy)

center = centroid(focus1, focus2)
linear_ecc = point_point_distance(focus1, focus2) / 2
major_axis_dir = (focus2[0] - focus1[0], focus2[1] - focus1[1])

print("  Semi-major axis length:")
print("  r₁ = (dist(f₁, p) + dist(f₂, point)) / 2\n")

r1 = symbols("r1")
r2 = sqrt(r1**2 - linear_ecc**2)

ellipse = ellipse(center, r1, r2, r1_direction=major_axis_dir)
ellipse /= gcd(list(ellipse))
ellipse = ellipse.applyfunc(cancel)

print("  Conic matrix:\n")
println_indented(ellipse)

print("This is equivalent to the matrix computed by ConicFromFociAndRadius:\n")

equiv_ellipse = expand(conic_from_foci_and_radius(focus1, focus2, r1))
println_indented(equiv_ellipse)

assert ellipse == equiv_ellipse
