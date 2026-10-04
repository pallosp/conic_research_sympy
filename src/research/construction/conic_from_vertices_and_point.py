#!/usr/bin/env python

"""Computes the conic whose vertices are `center ± v` and which contains the
point `p`.

In the results `r = |v|` is the distance between the center and the vertices.
"""

from sympy import (
    Eq,
    MatAdd,
    MatMul,
    Matrix,
    MatrixSymbol,
    Symbol,
    cancel,
    expand,
    factor,
    fraction,
    solve,
    sqrt,
    symbols,
    zeros,
)

from lib.central_conic import conic_from_center_and_points
from lib.circle import circle
from lib.degenerate_conic import line_pair_conic
from lib.line import line_between, line_through_point
from lib.point import ORIGIN, point_to_xy
from lib.transform import reflect_to_line, transform_point
from research.sympy_utils import println_indented

vx, vy, cx, cy, x, y = symbols("vx vy cx cy x y")
r = Symbol("r", positive=True)  # |v|
v = Matrix([vx, vy])
p = (x, y)
radius = sqrt(vx**2 + vy**2)


def conic_formula(center: tuple | Matrix) -> MatAdd:
    """Derives the conic with vertices at `center ± v` through `p` as an
    unevaluated sum of two matrices.

    The result is expressed with `r = |v|`: `vx² + vy²` is replaced with `r²`.
    """
    vertex1 = point_to_xy(center) + v
    vertex2 = point_to_xy(center) - v

    # The line through the vertices is an axis of the conic, therefore the
    # reflection of `p` to it is also on the conic.
    p_reflected = transform_point(p, reflect_to_line(line_between(vertex1, vertex2)))
    conic = conic_from_center_and_points(center, vertex1, p, p_reflected)

    # The conic is a member of the pencil spanned by
    #  - the pair of tangents at the vertices, and
    #  - the circle through the vertices.
    tangents = line_pair_conic(
        line_through_point(vertex1, normal=v),
        line_through_point(vertex2, normal=v),
    ).applyfunc(factor)
    vertex_circle = circle(center, radius)

    a, b = symbols("a b")
    equations = (conic - a * tangents - b * vertex_circle).applyfunc(expand)
    solution = solve(list(equations), [a, b])
    # Only the ratio of the coefficients matters
    circle_coeff, tangents_coeff = fraction(cancel(solution[b] / solution[a]))

    formula = MatAdd(
        MatMul(tangents, factor(tangents_coeff)),
        MatMul(vertex_circle, factor(circle_coeff)),
    ).subs(vx**2 + vy**2, r**2)

    # The formula is proportional to the conic computed from the center and
    # three points
    conic_from_formula = Matrix(formula.doit()).subs(r, radius)
    scale = cancel(conic[0] / conic_from_formula[0])
    assert (conic - scale * conic_from_formula).applyfunc(cancel) == zeros(3, 3)

    return formula


print("\nConic with vertices at center ± v through a point (r = |v|):\n")

c = MatrixSymbol("C", 3, 3)

print("\nCentered at the origin:\n")
println_indented(Eq(c, conic_formula(ORIGIN)))

print("\nCentered at (cx, cy):\n")
println_indented(Eq(c, conic_formula((cx, cy))))
