#!/usr/bin/env python

from itertools import product

from sympy import (
    Abs,
    Expr,
    I,
    Matrix,
    Ne,
    Symbol,
    cos,
    exp,
    expand_complex,
    factor,
    gcd,
    log,
    sin,
    sqrt,
    symbols,
)
from sympy.abc import x, y

from lib.conic import IdealPoints, conic_from_poly
from lib.conic_direction import ConicNormFactor
from lib.hyperbola import asymptote_focal_axis_angle
from lib.intersection import conic_x_line
from lib.line import IDEAL_LINE
from lib.matrix import conic_matrix, is_nonzero_multiple
from lib.transform import rotate, transform_point
from research.sympy_utils import eq_chain, println_indented

HORIZONTAL_LINE = "-" * 80


def normalize(point: Matrix) -> Matrix:
    """Divides out the common factor of the coordinates."""
    return (point / gcd(list(point))).applyfunc(factor)


################################################################################
# Intersection with the ideal line
################################################################################

print()
print("Ideal points on a conic:")
print()

println_indented(conic_x_line(conic_matrix(*symbols("a,b,c,d,e,f")), IDEAL_LINE))

print("Potential NonzeroCross expansions based on the coefficients' signs:")
print()

shown = set()
for a_hint, b_hint, c_hint in product(
    [{"zero": True}, {"nonzero": True}],
    [{"negative": True}, {"zero": True}, {"positive": True}],
    [{"zero": True}, {"nonzero": True}],
):
    a = Symbol("a", **a_hint)
    b = Symbol("b", **b_hint)
    c = Symbol("c", **c_hint)
    if a.equals(0) and b.equals(0) and c.equals(0):
        continue
    ideal_points = conic_x_line(conic_matrix(a, b, c, *symbols("d,e,f")), IDEAL_LINE)
    if str(ideal_points) not in shown:
        shown.add(str(ideal_points))
        println_indented(ideal_points)

################################################################################
# Formula based on focal axis angle and radii
################################################################################

print(HORIZONTAL_LINE)
print()
print("Radius and focal axis angle based formula:")
print()

r1, r2, lin_ecc, axis_angle = symbols("r1 r2 l alpha")

# Asymptote-axis angle; see ./asymptote_angle.py
cos_aa_angle = r1 / lin_ecc
sin_aa_angle = I * r2 / lin_ecc

# First ideal point angle: axis_angle + aa_angle
ideal_point_1 = Matrix(
    [
        cos(axis_angle) * cos_aa_angle - sin(axis_angle) * sin_aa_angle,
        cos(axis_angle) * sin_aa_angle + sin(axis_angle) * cos_aa_angle,
        0,
    ]
)

# Second ideal point angle: axis_angle - aa_angle
ideal_point_2 = ideal_point_1.subs(sin_aa_angle, -sin_aa_angle)

println_indented((normalize(ideal_point_1), normalize(ideal_point_2)))

################################################################################
# Branch-free, asymptote direction based formula
################################################################################

print(HORIZONTAL_LINE)
print()
print("Branch-free, asymptote direction based formula with complex coordinates:")
print()

a, c, d, e, f = symbols("a c d e f", real=True)
b = Symbol("b", real=True)
conic = conic_matrix(a, b, c, d, e, f)

eigen_minus, eigen_plus = symbols("lambda^- lambda^+", real=True)


def ideal_points_from_asymptotes(conic: Matrix) -> tuple[Matrix, Matrix]:
    """Rotates the focal axis direction by the ± asymptote-axis angle.

    The intermediate expressions are written with the eigenvalues `eigen_plus`
    and `eigen_minus` of the conic's quadratic part. They cancel out in the
    result.
    """
    a, _, _, b, c, _, _, _, _ = conic
    eigen_diff_square = (a - c) ** 2 + 4 * b**2
    norm = ConicNormFactor(conic)

    # Same as `focal_axis_direction(conic)`, but expressed with trigonometric
    # functions (`cos(atan2(2b, a-c)/2)`) instead of `Piecewise`, which breaks
    # `gcd` and `simplify` ("Piecewise generators do not make sense").
    axis_dir = Matrix([*sqrt(norm * (a - c + 2 * I * b)).simplify().as_real_imag(), 0])
    angle = asymptote_focal_axis_angle(conic)

    def eliminate_trig(coord: Expr) -> Expr:
        return (
            coord.subs(norm, (eigen_plus - eigen_minus) / Abs(eigen_plus - eigen_minus))
            .subs(eigen_diff_square.expand(), (eigen_plus - eigen_minus) ** 2)
            .subs(a + c, eigen_plus + eigen_minus)
            .rewrite(log)
            .factor(deep=True)
            .subs(eigen_diff_square.expand(), (eigen_plus - eigen_minus) ** 2)
            .rewrite(exp)
            .simplify()
        )

    ret = []
    for rotation in rotate(angle), rotate(-angle):
        point = transform_point(axis_dir, rotation).applyfunc(eliminate_trig)
        point = normalize(point).subs(1 / (eigen_plus - eigen_minus), 1)

        # Factor the coordinates as polynomials of sqrt(λ⁺) and sqrt(-λ⁻)
        sqrt_eigen_plus, sqrt_minus_eigen_minus = symbols("sep smem")
        point = (
            point.subs(sqrt(eigen_plus), sqrt_eigen_plus)
            .subs(sqrt(-eigen_minus), sqrt_minus_eigen_minus)
            .subs(eigen_plus, sqrt_eigen_plus**2)
            .subs(eigen_minus, -(sqrt_minus_eigen_minus**2))
        )
        point = normalize(point)

        # Substitute back, and use that λ⁺ + λ⁻ = a + c and λ⁺λ⁻ = ac - b²
        point = (
            point.subs(sqrt_eigen_plus, sqrt(eigen_plus))
            .subs(sqrt_minus_eigen_minus, sqrt(-eigen_minus))
            .subs(eigen_plus + eigen_minus, a + c)
            .subs(sqrt(eigen_plus) * sqrt(-eigen_minus), I * sqrt(a * c - b * b))
            .expand()
        )
        ret.append(point / 2)

    return ret[0], ret[1]


ip_formulae = ideal_points_from_asymptotes(conic)

println_indented(ip_formulae)

print("Verification for concrete conics:")
print()

conics = [
    conic_from_poly(x * x - y * y - 1),
    conic_from_poly(x * x - y * y + 1),
    conic_from_poly(x * y - 1),
    conic_from_poly(x * y + 1),
    conic_from_poly(2 * x * x + y * y - 1),
    conic_from_poly(2 * x * x + y * y + 1),
    conic_from_poly(x * x - y),
    conic_from_poly(x * x + 2 * y * y - 1),
    conic_from_poly(x * x - x * y - 1),
]

for conic_example in conics:
    ip1 = IdealPoints(conic_example)
    ip2 = tuple(
        expand_complex(formula.subs(zip(conic, conic_example, strict=True)))
        for formula in ip_formulae
    )
    assert any(is_nonzero_multiple(ip1[0], p) for p in ip2)
    assert any(is_nonzero_multiple(ip1[1], p) for p in ip2)
    println_indented(eq_chain(ip1, ip2))

print("The formula breaks down for circles:")
print()

circle = conic_from_poly(x * x + y * y - 1)
expected = IdealPoints(circle)[0]
actual = ip_formulae[0].subs(zip(conic, circle, strict=True))
println_indented(Ne(expected, actual, evaluate=False))
