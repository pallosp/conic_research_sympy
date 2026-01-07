#!/usr/bin/env python

from sympy import Abs, Matrix, Ne, exp, expand_complex, gcd, log, sqrt, symbols
from sympy.abc import x, y

from lib.conic import IdealPoints, conic_from_poly
from lib.conic_direction import ConicNormFactor, focal_axis_direction
from lib.hyperbola import asymptote_focal_axis_angle
from lib.intersection import conic_x_line
from lib.line import IDEAL_LINE
from lib.matrix import conic_matrix, is_nonzero_multiple
from lib.transform import rotate, transform_point
from research.sympy_utils import eq_chain, println_indented

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

    eigen_minus, eigen_plus = symbols("lambda^- lambda^+", real=True)
    a, _, _, b, c, _, _, _, _ = conic
    eigen_diff_square = (a - c) ** 2 + 4 * b**2

    def simplify_coord(point: Matrix) -> Matrix:
        eigen_diff = symbols("eigen_diff", nonzero=True)
        return point.applyfunc(
            lambda coord: (
                coord.subs(
                    ConicNormFactor(conic),
                    (eigen_plus - eigen_minus) / Abs(eigen_plus - eigen_minus),
                )
                .subs(eigen_diff_square.expand(), (eigen_plus - eigen_minus) ** 2)
                .subs(eigen_diff_square, (eigen_plus - eigen_minus) ** 2)
                .subs(a + c, eigen_plus + eigen_minus)
                .factor(deep=True)
                .rewrite(log)
                .factor(deep=True)
                .subs(eigen_diff_square.expand(), (eigen_plus - eigen_minus) ** 2)
                .subs(eigen_plus - eigen_minus, eigen_diff)
                .rewrite(exp)
                .simplify()
                .subs(eigen_diff, eigen_plus - eigen_minus)
            ),
        )

    ret = []
    for point in ideal_point_1, ideal_point_2:
        pt = simplify_coord(point)
        pt /= gcd(*pt[:2]).factor()
        pt = pt.subs(1 / (eigen_plus - eigen_minus), 1)
        if expand:
            pt = pt.subs(eigen_plus, (a + c + sqrt(eigen_diff_square)) / 2).subs(
                eigen_minus,
                (a + c - sqrt(eigen_diff_square)) / 2,
            )

        ret.append(pt)

    return tuple(ret)


conic = conic_matrix(*symbols("a b c d e f", real=True))
println_indented(ideal_points_from_asymptotes(conic, expand=False)[0])

print("Verification for concrete conics:\n")

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

print("The formula breaks down for circles:\n")

circle = conic_from_poly(x * x + y * y - 1)
expected = IdealPoints(circle)[0]
actual = ip_formula[0].subs(zip(conic, circle, strict=True))
println_indented(Ne(expected, actual, evaluate=False))
