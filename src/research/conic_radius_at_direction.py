#!/usr/bin/env python

from sympy import (
    Determinant,
    Eq,
    cos,
    oo,
    sin,
    solve,
    sqrt,
    symbols,
)
from sympy.simplify.fu import TR8

from lib.central_conic import conic_center
from lib.distance import point_point_distance
from lib.intersection import conic_x_line
from lib.line import line_through_point
from lib.matrix import conic_matrix
from research.sympy_utils import println_indented

HORIZONTAL_LINE = "=" * 80

################################################################################

print("\nLength of the conic radius at an angle of θ with horizontal:\n")

a, b, c, d, e, f = symbols("a b c d e f", real=True)
conic = conic_matrix(a, b, c, d, e, f)
center = conic_center(conic)
radius_angle = symbols("theta", real=True)
radius_dir = (cos(radius_angle), sin(radius_angle))
line_through_center = line_through_point(center, direction=radius_dir)
intersections = conic_x_line(conic, line_through_center)
intersection = intersections.args[0].row(2)

radius = (
    point_point_distance(center, intersection)
    .factor(deep=True)
    .trigsimp()
    .subs(conic.det(), Determinant(conic))
)

# Normalize the radius to be positive or a positive multiple of 𝑖
radius = sqrt(radius**2)

symbolic_radius_function = symbols("r(θ)")

println_indented(Eq(symbolic_radius_function, radius))

################################################################################

print(HORIZONTAL_LINE)

print("\nExpressed with the trigonometric functions of 2⋅θ:\n")

radius = TR8(radius).factor(deep=True).collect(cos(2 * radius_angle))
println_indented(Eq(symbolic_radius_function, radius))

################################################################################

print(HORIZONTAL_LINE)

print("\nLonger radius (in terms of absolute value):\n")

angle_subexpr = 2 * b * sin(2 * radius_angle) + (a - c) * cos(2 * radius_angle)
println_indented(radius.subs(angle_subexpr, -sqrt((2 * b) ** 2 + (a - c) ** 2)))

################################################################################

print(HORIZONTAL_LINE)

print("\nShorter radius (in terms of absolute value):\n")

println_indented(radius.subs(angle_subexpr, sqrt((2 * b) ** 2 + (a - c) ** 2)))

################################################################################

print(HORIZONTAL_LINE)

print("\nAngle of the asymptotes:\n")

println_indented(Eq(symbolic_radius_function, oo))

factor_with_angle = (
    a + c + (a - c) * cos(2 * radius_angle) + 2 * b * sin(2 * radius_angle)
)
println_indented(Eq(factor_with_angle, 0))

solutions = solve(factor_with_angle, radius_angle)
println_indented(Eq(symbols("theta_1"), solutions[0]))
println_indented(Eq(symbols("theta_2"), solutions[1]))

################################################################################
