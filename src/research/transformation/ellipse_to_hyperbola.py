#!/usr/bin/env python

from sympy import I, Matrix, factor, gcd, pi, sqrt, symbols

from lib.central_conic import (
    center_to_covertex_vector,
    center_to_vertex_vector,
    conic_center,
)
from lib.distance import point_point_distance
from lib.ellipse import ellipse
from lib.point import centroid, point_to_xy
from lib.transform import (
    homography_from_samples,
    rotate,
    scale,
    transform_conic,
    transform_point,
)
from research.sympy_utils import println_indented

print(
    "\nHomography that maps an ellipse to a hyperbola, preserving the (xᵢ, yᵢ) "
    "vertices and the absolute value of the secondary radius:\n"
)

vertex_1 = Matrix(symbols("x1 y1", real=True))
vertex_2 = Matrix(symbols("x2 y2", real=True))
r2 = symbols("r2", positive=True)

r1 = point_point_distance(vertex_1, vertex_2) / 2
center = centroid(vertex_1, vertex_2)
covertex_1 = transform_point(vertex_1, rotate(pi / 2, center) * scale(r2 / r1, center))
covertex_2 = transform_point(vertex_2, rotate(pi / 2, center) * scale(r2 / r1, center))
lin_ecc = sqrt(r1**2 - r2**2)

# Angle between the focal axis and the x-axis

cos_a, sin_a = (vertex_1 - center) / r1

# Angle between the focal axis and the asymptotes

cos_b, sin_b = r1 / lin_ecc, r2 / lin_ecc

# Angle between the ideal points and the x-axis

cos_a_plus_b = cos_a * cos_b - sin_a * sin_b
sin_a_plus_b = sin_a * cos_b + cos_a * sin_b

cos_a_minus_b = cos_a * cos_b + sin_a * sin_b
sin_a_minus_b = sin_a * cos_b - cos_a * sin_b

ideal_point_1 = Matrix([cos_a_plus_b, sin_a_plus_b, 0])
ideal_point_2 = Matrix([cos_a_minus_b, sin_a_minus_b, 0])

# Ellipse -> hyperbola transformation matrix

t = homography_from_samples(
    [vertex_1, vertex_2, covertex_1, covertex_2],
    [vertex_1, vertex_2, ideal_point_1, ideal_point_2],
)
t = t.applyfunc(lambda el: el.factor(deep=True))
t = t / gcd(list(t)).factor()

println_indented(t)

print("Expressed with the center point and the center-to-vertex vector:\n")

cx, cy, vx, vy = symbols("cx cy vx vy")
println_indented(
    (t / 4)
    .subs(vertex_1[0], cx + vx)
    .subs(vertex_1[1], cy + vy)
    .subs(vertex_2[0], cx - vx)
    .subs(vertex_2[1], cy - vy)
    .applyfunc(factor)
)

print("Interesting fact: the transformation is independent on the secondary radius.\n")

# Verify correctness

h_ellipse = ellipse((4, 3), 2, 1)
hyperbola = transform_conic(
    h_ellipse,
    t.subs(zip(vertex_1, (6, 3), strict=True)).subs(zip(vertex_2, (2, 3), strict=True)),
)

assert point_to_xy(conic_center(hyperbola)) == Matrix([4, 3])
assert center_to_vertex_vector(hyperbola) == Matrix([2, 0])
assert center_to_covertex_vector(hyperbola) == Matrix([0, I])
