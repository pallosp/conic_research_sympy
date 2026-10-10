#!/usr/bin/env python

from sympy import expand, factor, symbols

from lib.matrix import conic_matrix
from research.sympy_utils import println_indented

# Circle from center point c=(cx, cy) and a finite point p=(px, py)
# Equation: (x-cx)² + (y-cy)² = r², where r² = (px-cx)² + (py-cy)²

cx, cy, px, py = symbols("cx cy px py")
r_squared = (px - cx) ** 2 + (py - cy) ** 2
# The sign is chosen so that the conic's value at the center is r² >= 0.
circle = conic_matrix(-1, 0, -1, cx, cy, r_squared - cx**2 - cy**2)

print("\nCircle from center point and a finite point on it:\n")
println_indented(circle.applyfunc(expand))

# Extend the formula to projective points p=(px, py, pz) by substituting
# px/pz and py/pz. Multiplying the matrix by pz² clears the denominators.

pz = symbols("pz")
circle = (circle.subs({px: px / pz, py: py / pz}) * pz**2).applyfunc(
    lambda e: factor(expand(e))
)

print("\nCircle from center point and a projective point p on it:\n")
println_indented(circle)

# When p is ideal, the radius is infinite, and the circle degenerates into the
# double ideal line (unless px² + py² = 0), independently of the center.

print("\nWhen p is ideal:\n")
ideal_p_circle = circle.subs({pz: 0})
println_indented(ideal_p_circle)
