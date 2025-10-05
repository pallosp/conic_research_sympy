#!/usr/bin/env python

from sympy import cancel, gcd, sqrt, symbols

from lib.distance import PointPointDistance
from lib.ellipse import Ellipse
from lib.point import Centroid
from research.util import println_indented

print("\nEllipse from its focus points and an incident point:\n")

cx, cy = symbols("cx cy")  # center
dx, dy = symbols("dx dy")  # center -> focus vector

focus1 = (cx + dx, cy + dy)
focus2 = (cx - dx, cy - dy)

center = Centroid(focus1, focus2)
linear_ecc = PointPointDistance(focus1, focus2) / 2
major_axis_dir = (focus2[0] - focus1[0], focus2[1] - focus1[1])

print("  Semi-major axis length:")
print("  r₁ = (dist(f₁, p) + dist(f₂, point)) / 2\n")

r1 = symbols("r1")
r2 = sqrt(r1**2 - linear_ecc**2)

ellipse = Ellipse(center, r1, r2, r1_direction=major_axis_dir)
ellipse /= gcd(list(ellipse))
ellipse = ellipse.applyfunc(cancel)

print("  Conic matrix:\n")
println_indented(ellipse)
