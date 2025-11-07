#!/usr/bin/env python

from sympy import (
    Abs,
    Matrix,
    exp,
    expand_complex,
    gcd,
    log,
    sqrt,
    symbols,
)
from sympy.abc import x, y

from lib.conic import IdealPoints, conic_from_poly, focal_axis_direction
from lib.conic_classes import ConicNormFactor
from lib.hyperbola import asymptote_focal_axis_angle
from lib.intersection import conic_x_line
from lib.line import IDEAL_LINE
from lib.matrix import conic_matrix, is_nonzero_multiple
from lib.sympy_utils import eq_chain
from lib.transform import rotate, transform_point
from research.util import println_indented

print("\nIdeal points on a conic:\n")

println_indented(conic_x_line(conic_matrix(*symbols("a,b,c,d,e,f")), IDEAL_LINE))

print("\nPotential formulae based on the coefficients' signs:\n")

solutions = set()
for a_assumption in [{"zero": True}, {"nonzero": True}]:
    for b_assumption in [{"negative": True}, {"zero": True}, {"positive": True}]:
        for c_assumption in [{"zero": True}, {"nonzero": True}]:
            a = symbols("a", **a_assumption)
            b = symbols("b", **b_assumption)
            c = symbols("c", **c_assumption)
            if a.equals(0) and b.equals(0) and c.equals(0):
                continue
            d, e, f = symbols("d,e,f")
            conic = conic_matrix(a, b, c, d, e, f)
            ideal_points = conic_x_line(conic, IDEAL_LINE)
            if str(ideal_points) not in solutions:
                solutions.add(str(ideal_points))
                println_indented(ideal_points)

print("Asymptote direction based formula:\n")


def ideal_points_from_asymptotes(
    conic: Matrix,
    *,
    expand: bool = False,
) -> tuple[Matrix, Matrix]:
    axis_dir = focal_axis_direction(conic)
    rotation1 = rotate(asymptote_focal_axis_angle(conic))
    rotation2 = rotate(-asymptote_focal_axis_angle(conic))
    ideal_point_1 = transform_point(axis_dir, rotation1)
    ideal_point_2 = transform_point(axis_dir, rotation2)

    n, p = symbols("lambda^- lambda^+", real=True)
    p_minus_n = symbols("pmn", nonzero=True)
    a, _, _, b, c, _, _, _, _ = conic
    eigen_diff_square = (a - c) ** 2 + 4 * b**2

    def simplify_coord(point: Matrix) -> Matrix:
        return point.applyfunc(
            lambda coord: (
                coord.subs(ConicNormFactor(conic), (p - n) / Abs(p - n))
                .subs(eigen_diff_square.expand(), (p - n) ** 2)
                .subs(eigen_diff_square, (p - n) ** 2)
                .subs(a + c, p + n)
                .factor(deep=True)
                .rewrite(log)
                .factor(deep=True)
                .subs(eigen_diff_square.expand(), (p - n) ** 2)
                .subs(p - n, p_minus_n)
                .rewrite(exp)
                .simplify()
            ),
        )

    ret = []
    for point in ideal_point_1, ideal_point_2:
        pt = simplify_coord(point)
        pt /= gcd(*pt[:2]).factor() / 2 / p_minus_n
        if expand:
            pt = (
                pt.subs(p_minus_n, p - n)
                .subs(p, (a + c + sqrt(eigen_diff_square)) / 2)
                .subs(n, (a + c - sqrt(eigen_diff_square)) / 2)
            )

        ret.append(pt)

    return tuple(ret)


conic = conic_matrix(*symbols("a b c d e f", real=True))
println_indented(ideal_points_from_asymptotes(conic, expand=False)[0])

print("\nVerification for concrete conics:\n")

conics = [
    conic_from_poly(x * x - y * y - 1),
    conic_from_poly(x * x - y * y + 1),
    conic_from_poly(x * y - 1),
    conic_from_poly(2 * x * x + y * y - 1),
    conic_from_poly(2 * x * x + y * y + 1),
    conic_from_poly(x * x - y),
]
ip_formula = ideal_points_from_asymptotes(conic, expand=True)
for conic_example in conics:
    ip1 = IdealPoints(conic_example)
    ip2 = tuple(
        expand_complex(i.subs(zip(conic, conic_example, strict=True)))
        for i in ip_formula
    )
    assert any(is_nonzero_multiple(ip1[0], p) for p in ip2)
    assert any(is_nonzero_multiple(ip1[1], p) for p in ip2)
    println_indented(eq_chain(ip1, ip2))
