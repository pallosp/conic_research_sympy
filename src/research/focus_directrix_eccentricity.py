#!/usr/bin/env python

from textwrap import indent

from sympy import (
    Eq,
    Matrix,
    Piecewise,
    Pow,
    Q,
    atan2,
    cos,
    factor,
    pretty,
    sign,
    simplify,
    sin,
    solve,
    sqrt,
    symbols,
)

from lib.conic import ConicFromFocusAndDirectrix
from lib.matrix import ConicMatrix
from lib.sympy_utils import AddEq, DivEq, FactorRadicals, MulEq, SubEq, SwapEq


def print_indented(expr: object) -> None:
    print(indent(pretty(expr), "  "))


def println_indented(expr: object) -> None:
    print(indent(pretty(expr), "  ") + "\n")


fx, fy, a, b, c, ecc, L = symbols("x,y,a,b,c,e,lambda", real=True)

focus = (fx, fy)
directrix = Matrix([a, b, c])
conic = ConicFromFocusAndDirectrix(focus, directrix, ecc) * L


print("\nEquations to determine the foci, directrices and eccentricity:\n")

A, B, C, D, E, F = symbols("A,B,C,D,E,F", real=True)
coeff_eq = [
    Eq(A, conic[0]),
    Eq(B, conic[3]),
    Eq(C, conic[4]),
    Eq(D, conic[6]),
    Eq(E, conic[7]),
    Eq(F, conic[8]),
]
for eq in coeff_eq:
    print_indented(eq)

print("\nSince a² + b² = 1:\n")

for i in range(len(coeff_eq)):
    coeff_eq[i] = coeff_eq[i].subs(a * a + b * b, 1)
    print_indented(coeff_eq[i])

print("\nDeterminant of the conic matrix:\n")

det = symbols("det", nonzero=True)
det_eq = Eq(det, conic.det().factor().subs(a * a + b * b, 1))
print_indented(det_eq)

print("\nThis implies that sign(λ) = sign(det).")
print("Note that the formula breaks down for circles.\n")

print("Calculate the eigenvalues of the top left submatrix.\n")

min_eigen_abc, max_eigen_abc = Matrix([[A, B], [B, C]]).eigenvals()
min_eigen_fde, max_eigen_fde = (
    e.subs([(eq.lhs, eq.rhs) for eq in coeff_eq])
    .factor()
    .replace(Pow, lambda base, exp: Pow(base.factor(), exp))
    .subs(b * b, 1 - a * a)
    .factor()
    .collect(ecc)
    for e in [min_eigen_abc, max_eigen_abc]
)


min_eigen, max_eigen = symbols("min_eigen,max_eigen", real=True, nonzero=True)
eigenvalue_min_eq_abc = Eq(min_eigen, min_eigen_abc)
eigenvalue_max_eq_abc = Eq(max_eigen, max_eigen_abc)
eigenvalue_sum_eq_abc = AddEq(eigenvalue_min_eq_abc, eigenvalue_max_eq_abc)
eigenvalue_prod_eq_abc = simplify(MulEq(eigenvalue_max_eq_abc, eigenvalue_min_eq_abc))
eigenvalue_diff_eq_abc = SubEq(eigenvalue_max_eq_abc, eigenvalue_min_eq_abc)
eigenvalue_min_eq_fde = Eq(min_eigen, min_eigen_fde)
eigenvalue_max_eq_fde = Eq(max_eigen, max_eigen_fde)
eigenvalue_sum_eq_fde = factor(AddEq(eigenvalue_min_eq_fde, eigenvalue_max_eq_fde))
eigenvalue_diff_eq_fde = factor(SubEq(eigenvalue_max_eq_fde, eigenvalue_min_eq_fde))

eigenvalue_eqs = [
    eigenvalue_min_eq_abc,
    eigenvalue_max_eq_abc,
    eigenvalue_sum_eq_abc,
    eigenvalue_prod_eq_abc,
    eigenvalue_diff_eq_abc,
    eigenvalue_min_eq_fde,
    eigenvalue_max_eq_fde,
    eigenvalue_sum_eq_fde,
    eigenvalue_diff_eq_fde,
]

for eq in eigenvalue_eqs:
    println_indented(eq)

print("As sign(λ) = sign(det):\n")

eigenvalue_diff_eq_fde = Eq(
    (max_eigen - min_eigen) * sign(det, evaluate=False),
    ecc**2 * L,
)
println_indented(eigenvalue_diff_eq_fde)

lambda_eq = DivEq(SubEq(eigenvalue_diff_eq_fde, eigenvalue_sum_eq_fde), 2)
lambda_eq = SwapEq(factor(lambda_eq))

println_indented(lambda_eq)

lambda_eq_piecewise = Eq(L, Piecewise((-min_eigen, det > 0), (-max_eigen, True)))
println_indented(lambda_eq_piecewise)

eigen_to_abc = {min_eigen: min_eigen_abc, max_eigen: max_eigen_abc}
lambda_eq_abc = factor(lambda_eq.subs(eigen_to_abc))
println_indented(lambda_eq_abc)

ecc_square_eq = SwapEq(DivEq(eigenvalue_diff_eq_fde, L))
println_indented(ecc_square_eq)

ecc_square_eq_abc = ecc_square_eq.subs(L, lambda_eq_abc.rhs).subs(eigen_to_abc)
println_indented(ecc_square_eq_abc)

print("Directrix normal vector a.k.a. focal axis direction:\n")

aa_eq = Eq(a**2, solve(coeff_eq[0], a**2)[0])
bb_eq = Eq(b**2, solve(coeff_eq[2], b**2)[0])
ab_eq = SwapEq(DivEq(coeff_eq[1], ecc**2 * L))
println_indented(aa_eq)
println_indented(bb_eq)
println_indented(ab_eq)

print("As a² + b² = 1, substitute a = cos(θ), b = sin(θ)\n")

theta = symbols("theta")
println_indented(SubEq(aa_eq, bb_eq))
println_indented(SubEq(aa_eq, bb_eq).subs({a: cos(theta), b: sin(theta)}))
cos2_eq = SubEq(aa_eq, bb_eq).subs({a: cos(theta), b: sin(theta)}).simplify()
println_indented(cos2_eq)
println_indented(MulEq(ab_eq, 2))
println_indented(MulEq(ab_eq, 2).subs({a: cos(theta), b: sin(theta)}))
sin2_eq = MulEq(ab_eq, 2).subs({a: cos(theta), b: sin(theta)}).simplify()
println_indented(sin2_eq)
theta_eq = Eq(theta, atan2(sin2_eq.rhs, cos2_eq.rhs) / 2)
println_indented(theta_eq)
println_indented(theta_eq.subs(ecc**2 * L, sign(ecc**2 * L)))
println_indented(theta_eq.subs(ecc**2 * L, 1 / sign(det)))

print("\nDirectrix equation:\n")

println_indented(coeff_eq[3])
fx_formula = solve(coeff_eq[3], fx)[0]
println_indented(Eq(fx, fx_formula))

println_indented(coeff_eq[4])
fy_formula = solve(coeff_eq[4], fy)[0]
println_indented(Eq(fy, fy_formula))

f_eq = coeff_eq[5]
println_indented(f_eq)
f_eq = f_eq.subs({fx: fx_formula, fy: fy_formula})
println_indented(f_eq)

print("\nSolutions when there are two directrices:\n")

directrix_c_values = solve(f_eq, c)
assert len(directrix_c_values) == 2

for det_assumption in Q.positive, Q.negative:
    print(f"\n  when det{'>0' if det_assumption == Q.positive else '<0'}\n")
    for index, c in enumerate(directrix_c_values):
        c_simplified = (
            c.subs(a * a, aa_eq.rhs)
            .subs(b * b, bb_eq.rhs)
            .subs(a * b, ab_eq.rhs)
            .subs(ecc, sqrt(ecc_square_eq.rhs))
            .subs(L, lambda_eq.rhs)
            .simplify()
            .refine(det_assumption(det))
            .factor()
            .collect(F)
            .subs(A + C, min_eigen + max_eigen)
            .simplify()
        )
        c_simplified = (
            FactorRadicals(c_simplified)
            .subs(max_eigen * min_eigen, eigenvalue_prod_eq_abc.rhs)
            .subs(max_eigen - min_eigen, eigenvalue_diff_eq_abc.rhs)
        )
        c_simplified = (
            FactorRadicals(c_simplified)
            .subs(ConicMatrix(A, B, C, D, E, F).det(), det)
            .subs(((A - C) ** 2).expand(), (A - C) ** 2)
            .refine(det_assumption(det))
        )
        println_indented(Eq(symbols(f"c{index+1}"), c_simplified))

print("\nSpecial case: circle\n")

print("  ∜(4B²+(A-C)²) = 0 ==> the c coefficent of the directrix is not valid.\n")
print("  Multiplying all coefficients with that expression however will yield ")
print("  a'=b'=0, c'≠0, which represents the ideal line.\n")
