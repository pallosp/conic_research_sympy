#!/usr/bin/env python

from sympy import (
    Abs,
    Eq,
    MatMul,
    Matrix,
    cos,
    factor,
    pi,
    sign,
    simplify,
    sin,
    solve,
    symbols,
)

from lib.central_conic import (
    conic_center,
    radius_in_direction,
)
from lib.ellipse import ellipse
from lib.matrix import conic_matrix
from lib.polar_conic import conic_from_polar_matrix, curvature_sign_at_angle
from lib.transform import homography_from_samples
from research.sympy_utils import eq_chain, println_indented

# ---------------------------------------------------------------------------
# First approach: start from the endpoint of the horizontal diameter
# ---------------------------------------------------------------------------

a, b, c, d, e, f = symbols("a b c d e f", real=True)
cx, cy = symbols("cx cy", real=True)
r_horizontal = symbols("rh", real=True)

print("Ellipse in polar form:")
print()

polar_ellipse = Matrix(3, 3, symbols("A B C D E F G H I", real=True))
println_indented(eq_chain(symbols("P"), polar_ellipse))

print(
    "The secants between α and π+α intersect at the point [C, F, I]ᵀ.\n"
    "In a symmetric polar representation of an ellipse, this intersection\n"
    "must coincide with the ellipse center."
)
print()

polar_ellipse = polar_ellipse.subs(zip(polar_ellipse.col(2), (cx, cy, 1), strict=True))
println_indented(eq_chain(symbols("P"), polar_ellipse))

print(
    "The horizontal line through the center intersects the ellipse in two points.\n"
    "These correspond to the angles 0 and π on the polar curve."
)
print()

print(
    "The secants between α and −α intersect at [B, E, H]ᵀ.\n"
    "By symmetry, this intersection must also be an ideal point."
)
print()

A, B, _, D, E, _, G, H, _ = polar_ellipse
polar_ellipse = polar_ellipse.subs(H, 0)
println_indented(eq_chain(symbols("P"), polar_ellipse))

print(
    "The secants between α and π−α intersect at [A, D, G]ᵀ.\n"
    "The horizontal diameter (between 0 and π) is one such secant.\n"
    "By symmetry, [A, D, G]ᵀ must therefore be the ideal point in the (1, 0) direction."
)
print()

polar_ellipse_h = polar_ellipse.subs(D, 0).subs(G, 0)
println_indented(eq_chain(symbols("P"), polar_ellipse_h))

print(
    "Let the horizontal radius be rh.\n"
    "Then the x-coordinate of the point P(0) is cx + rh.\n"
    "\n"
    "P(0) = P·[cos 0, sin 0, 1]ᵀ = [A + cx, cy, 1]ᵀ = λ·[cx + rh, cy, 1]ᵀ\n"
    "This yields the solution λ = 1 and A = rh."
)
print()

polar_ellipse_h = polar_ellipse_h.subs(A, r_horizontal)
println_indented(eq_chain(symbols("P"), polar_ellipse_h))

print("The ellipse in Cartesian form:")
print()

cartesian_ellipse = conic_from_polar_matrix(polar_ellipse_h)
mu = symbols("mu", nonzero=True)

println_indented(
    eq_chain(
        cartesian_ellipse,
        MatMul(mu, conic_matrix(a, b, c, d, e, f)),
    )
)

print("Solving for the remaining parameters:")
print()

solutions = solve(
    [
        Eq(cartesian_ellipse[0], mu * a),
        Eq(cartesian_ellipse[1], mu * b),
        Eq(cartesian_ellipse[4], mu * c),
    ],
    (mu, B, E),
)

for solution in solutions:
    println_indented(eq_chain((mu, B, E), solution))

print("Final polar matrix (using the first solution):")
print()

polar_ellipse_h = polar_ellipse_h.subs(B, solutions[0][1]).subs(E, solutions[0][2])
println_indented(eq_chain(symbols("P"), polar_ellipse_h))

cartesian_conic = conic_matrix(a, b, c, d, e, f)
center = conic_center(cartesian_conic)
r_horizontal_value = radius_in_direction(cartesian_conic, angle=0)
discriminant = symbols("disc", positive=True)

polar_ellipse_h = (
    polar_ellipse_h.subs(zip((cx, cy), center, strict=True))
    .subs(r_horizontal, r_horizontal_value)
    .subs(cartesian_conic.det(), symbols("det"))
    .subs(a * c - b * b, discriminant)
    .subs(discriminant, a * c - b * b)
)

println_indented(eq_chain(symbols("P"), polar_ellipse_h))

print("Orientation (1 is counterclockwise, -1 is clockwise):")
print()

println_indented(curvature_sign_at_angle(polar_ellipse_h, 0).factor(deep=True))

# ---------------------------------------------------------------------------
# Second approach: start from the endpoint of the vertical diameter
# ---------------------------------------------------------------------------

print("-" * 80)
print()
print("Polar ellipse representation starting at the endpoint of the vertical diameter")
print()

polar_ellipse_v = polar_ellipse.subs(A, 0).subs(G, 0)
r_vertical = symbols("rv", positive=True)

# P(0) = P·[cos 0, sin 0, 1]ᵀ = [cx, D + cy, 1]ᵀ = λ·[cx, cy + rv, 1]ᵀ\n"
# Solution: λ=1, D=rv
polar_ellipse_v = polar_ellipse_v.subs(D, r_vertical)

cartesian_ellipse = conic_from_polar_matrix(polar_ellipse_v)
solutions = solve(
    [
        Eq(cartesian_ellipse[0], mu * a),
        Eq(cartesian_ellipse[1], mu * b),
        Eq(cartesian_ellipse[4], mu * c),
    ],
    (mu, B, E),
)
polar_ellipse_v = polar_ellipse_v.subs(B, solutions[0][1]).subs(E, solutions[0][2])

println_indented(eq_chain(symbols("P"), polar_ellipse_v))

cartesian_conic = conic_matrix(a, b, c, d, e, f)
center = conic_center(cartesian_conic)
r_vertical_value = radius_in_direction(cartesian_conic, angle=pi / 2)
discriminant = symbols("disc", positive=True)

polar_ellipse_v = (
    polar_ellipse_v.subs(zip((cx, cy), center, strict=True))
    .subs(r_vertical, r_vertical_value)
    .subs(cartesian_conic.det(), symbols("det"))
    .subs(a * c - b * b, discriminant)
    .subs(discriminant, a * c - b * b)
    .subs(Abs(b), b * sign(b))
)

println_indented(eq_chain(symbols("P"), polar_ellipse_v))

# ---------------------------------------------------------------------------
# Third approach: start from a vertex
# ---------------------------------------------------------------------------

print("-" * 80)
print()

print(
    "Alternative construction: map an ellipse vertex to angle = 0.\n"
    "Here, α denotes the angle between the focal axis and the horizontal."
)
print()

cx, cy, r1, r2, alpha = symbols("cx cy r1 r2 alpha")

t = homography_from_samples(
    [(1, 0), (0, 1), (-1, 0), (0, -1)],
    [
        (cx + r1 * cos(alpha), cy + r1 * sin(alpha)),
        (cx - r2 * sin(alpha), cy + r2 * cos(alpha)),
        (cx - r1 * cos(alpha), cy - r1 * sin(alpha)),
        (cx + r2 * sin(alpha), cy - r2 * cos(alpha)),
    ],
)

t /= t[8]
t = simplify(t)

c1 = conic_from_polar_matrix(t).applyfunc(factor).applyfunc(simplify)
c2 = ellipse((cx, cy), r1, r2, r1_angle=alpha).applyfunc(factor).applyfunc(simplify)

assert c1 == c2

println_indented(eq_chain(symbols("P"), simplify(t)))

print(
    "Although the resulting matrix appears simple,\n"
    "there is substantial algebraic complexity hidden in expressing\n"
    "cos(α) and sin(α) in terms of the Cartesian conic coefficients."
)
print()
