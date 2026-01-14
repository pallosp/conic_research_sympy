#!/usr/bin/env python3

from sympy import (
    Expr,
    Matrix,
    cos,
    expand_trig,
    factor,
    sin,
    symbols,
)

from lib.polar_conic import point_at_angle
from lib.transform import rotate
from research.sympy_utils import println_indented

a, b, c, d, e, f, g, h, i = symbols("a b c d e f g h i")
polar_conic = Matrix([[a, b, c], [d, e, f], [g, h, i]])

alpha, beta = symbols("alpha beta")
p0 = point_at_angle(polar_conic, alpha)
p1 = point_at_angle(polar_conic, alpha + beta).applyfunc(expand_trig)


def decode_polar_matrix_row(point_coord: Expr, angle: Expr) -> list[Expr]:
    """Converts a·cos(angle)+b·sin(angle)+c to [a, b, c]."""
    point_coord = point_coord.expand()
    components = point_coord.collect([cos(angle), sin(angle)], evaluate=False)
    return [components[cos(angle)], components[sin(angle)], components[1]]


polar_matrix_for_alpha = Matrix([decode_polar_matrix_row(c, alpha) for c in p1])
pre_rotation = (polar_conic.inv() * polar_matrix_for_alpha).expand().applyfunc(factor)

print("\nHow to rotate the points of a polar conic P along the curve?\n")

print("Matrix T, for which P·T·[cos α, sin α, 1]ᵀ = P·[cos α+β, sin α+β, 1]ᵀ:\n")
println_indented(pre_rotation)

assert p1.expand() == point_at_angle(polar_conic * pre_rotation, alpha).expand()
assert pre_rotation == rotate(beta)
