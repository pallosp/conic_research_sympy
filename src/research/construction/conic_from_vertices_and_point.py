#!/usr/bin/env python

"""Computes the conic whose vertices are `center ± v` and which contains the
point `p`.

In the results `r = |v|` is the distance between the center and the vertices.

The conics touching the vertex tangents at the two vertices form a pencil. The
results express the conic as a sum of two members of this pencil:
 - the pair of tangents at the vertices: the lines through `center + v` and
   `center - v` perpendicular to `v`, as a degenerate conic, and
 - either the circle whose diameter is the segment between the vertices, or
   the double line through the center in the direction of `v`, i.e. the line
   containing the vertices counted twice.
"""

from itertools import combinations

from sympy import (
    Add,
    Eq,
    Expr,
    MatAdd,
    MatMul,
    Matrix,
    MatrixSymbol,
    Mul,
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
from lib.degenerate_conic import double_line_conic, line_pair_conic
from lib.line import line_between, line_through_point
from lib.point import ORIGIN, point_to_xy
from lib.transform import reflect_to_line, transform_point
from research.sympy_utils import println_indented

vx, vy, cx, cy, x, y = symbols("vx vy cx cy x y")
r = Symbol("r", positive=True)  # |v|
v = Matrix([vx, vy])
p = (x, y)
radius = sqrt(vx**2 + vy**2)


def difference_of_squares(expr: Expr) -> Expr:
    """Rewrites `(a - b)(a + b)` subexpressions to `a² - b²`.

    `a` is the sum of the terms the two factors share.
    """

    def rewrite(product: Mul) -> Expr:
        factors = list(product.args)
        for i, j in combinations(range(len(factors)), 2):
            f, g = factors[i], factors[j]
            if not (isinstance(f, Add) and isinstance(g, Add)):
                continue
            common = set(f.args) & set(g.args)
            f_rest = f - Add(*common)
            g_rest = g - Add(*common)
            if common and f_rest != 0 and expand(f_rest + g_rest) == 0:
                factors[i] = Add(*common) ** 2 - g_rest**2
                del factors[j]
                return Mul(*factors)
        return product

    return expr.replace(lambda e: isinstance(e, Mul), rewrite)


def conic_formula(center: tuple | Matrix, *, double_line: bool = False) -> MatAdd:
    """Derives the conic with vertices at `center ± v` through `p` as an
    unevaluated sum of two matrices.

    The sum consists of two members of the pencil of conics touching the
    vertex tangents at the vertices:
     - the pair of tangents at the vertices: the lines through `center + v` and
       `center - v` perpendicular to `v`, as a degenerate conic, and
     - the circle centered at `center` with radius `|v|` (default), or the double
       line through `center` in the direction of `v` (`double_line=True`).

    The result is expressed with `r = |v|`: `vx² + vy²` is replaced with `r²`.
    """
    vertex1 = point_to_xy(center) + v
    vertex2 = point_to_xy(center) - v

    # The line through the vertices is an axis of the conic, therefore the
    # reflection of `p` to it is also on the conic.
    p_reflected = transform_point(p, reflect_to_line(line_between(vertex1, vertex2)))
    conic = conic_from_center_and_points(center, vertex1, p, p_reflected)

    # The conic touches the vertex tangents at the vertices. The conics with this
    # property form a pencil, so the conic is a combination of any two of
    #  - the pair of tangents at the vertices,
    #  - the circle with the vertices as the endpoints of a diameter,
    #  - the double line through the vertices.
    tangents = line_pair_conic(
        line_through_point(vertex1, normal=v),
        line_through_point(vertex2, normal=v),
    ).applyfunc(factor)
    if double_line:
        axis = line_through_point(center, direction=v)
        second_member = double_line_conic(axis).applyfunc(factor)
    else:
        second_member = -circle(center, radius)

    a, b = symbols("a b")
    equations = (conic - a * tangents - b * second_member).applyfunc(expand)
    solution = solve(list(equations), [a, b])
    # Only the ratio of the coefficients matters
    second_coeff, tangents_coeff = fraction(cancel(solution[b] / solution[a]))

    formula = MatAdd(
        MatMul(
            tangents.applyfunc(difference_of_squares),
            difference_of_squares(factor(tangents_coeff)),
        ),
        MatMul(
            second_member.applyfunc(difference_of_squares),
            difference_of_squares(factor(second_coeff)),
        ),
    ).subs(vx**2 + vy**2, r**2)

    # The formula is proportional to the conic computed from the center and
    # three points
    conic_from_formula = Matrix(formula.doit()).subs(r, radius)
    scale = cancel(conic[0] / conic_from_formula[0])
    assert (conic - scale * conic_from_formula).applyfunc(cancel) == zeros(3, 3)

    return formula


print("\nConic with vertices at center ± v through a point (r = |v|):\n")

c = MatrixSymbol("C", 3, 3)

for center, center_name in [(ORIGIN, "the origin"), ((cx, cy), "(cx, cy)")]:
    for double_line, member_name in [
        (False, "circle with the vertices as diameter endpoints"),
        (True, "double line through the vertices"),
    ]:
        print(
            f"\nCentered at {center_name}, pencil members: "
            f"tangents at the vertices and {member_name}:\n",
        )
        println_indented(Eq(c, conic_formula(center, double_line=double_line)))
