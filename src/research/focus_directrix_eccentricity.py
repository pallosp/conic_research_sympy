#!/usr/bin/env python

from textwrap import indent

from sympy import (
    Eq,
    Expr,
    Matrix,
    Piecewise,
    Pow,
    atan2,
    cos,
    factor,
    pretty,
    sign,
    sin,
    solve,
    symbols,
)

from lib.conic import ConicFromFocusAndDirectrix


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


def Add(*eqs: Eq) -> Eq:
    left = eqs[0].lhs
    right = eqs[0].rhs
    for term in eqs[1:]:
        left += term.lhs
        right += term.rhs
    return Eq(left, right)


def Sub(eq0: Eq, eq1: Eq) -> Eq:
    return Eq(eq0.lhs - eq1.lhs, eq0.rhs - eq1.rhs)


def Mul(eq: Eq, factor: Expr) -> Eq:
    return Eq(eq.lhs * factor, eq.rhs * factor)


def Div(eq: Eq, denom: Expr) -> Eq:
    return Eq(eq.lhs / denom, eq.rhs / denom)


def Swap(eq: Eq) -> Eq:
    return Eq(eq.rhs, eq.lhs)


min_eigen, max_eigen = symbols("min_eigen,max_eigen")
eigenvalue_min_eq = Eq(min_eigen, min_eigen_fde)
eigenvalue_max_eq = Eq(max_eigen, max_eigen_fde)
eigenvalue_sum_eq = factor(Add(eigenvalue_min_eq, eigenvalue_max_eq))
eigenvalue_diff_eq = factor(Sub(eigenvalue_max_eq, eigenvalue_min_eq))

eigenvalue_eqs = [
    eigenvalue_min_eq,
    eigenvalue_max_eq,
    eigenvalue_sum_eq,
    eigenvalue_diff_eq,
]

for eq in eigenvalue_eqs:
    println_indented(eq)

print("As sign(λ) = sign(det):\n")

eigenvalue_diff_eq = Eq((max_eigen - min_eigen) * sign(det, evaluate=False), ecc**2 * L)
println_indented(eigenvalue_diff_eq)

lambda_eq = Swap(factor(Div(Sub(eigenvalue_diff_eq, eigenvalue_sum_eq), 2)))
println_indented(lambda_eq)

lambda_eq_piecewise = Eq(L, Piecewise((-min_eigen, det > 0), (-max_eigen, True)))
println_indented(lambda_eq_piecewise)

eigen_to_abc = {min_eigen: min_eigen_abc, max_eigen: max_eigen_abc}
lambda_eq_abc = factor(lambda_eq.subs(eigen_to_abc))
println_indented(lambda_eq_abc)

ecc_square_eq = Swap(Div(eigenvalue_diff_eq, L))
println_indented(ecc_square_eq)

ecc_square_eq_abc = ecc_square_eq.subs(L, lambda_eq_abc.rhs).subs(eigen_to_abc)
println_indented(ecc_square_eq_abc)

print("Directrix normal vector a.k.a. focal axis direction:\n")

aa_eq = Eq(a**2, solve(coeff_eq[0], a**2)[0])
bb_eq = Eq(b**2, solve(coeff_eq[2], b**2)[0])
ab_eq = Swap(Div(coeff_eq[1], ecc**2 * L))
println_indented(aa_eq)
println_indented(bb_eq)
println_indented(ab_eq)

print("As a² + b² = 1, substitute a = cos(θ), b = sin(θ)\n")

theta = symbols("theta")
println_indented(Sub(aa_eq, bb_eq))
println_indented(Sub(aa_eq, bb_eq).subs({a: cos(theta), b: sin(theta)}))
cos2_eq = Sub(aa_eq, bb_eq).subs({a: cos(theta), b: sin(theta)}).simplify()
println_indented(cos2_eq)
println_indented(Mul(ab_eq, 2))
println_indented(Mul(ab_eq, 2).subs({a: cos(theta), b: sin(theta)}))
sin2_eq = Mul(ab_eq, 2).subs({a: cos(theta), b: sin(theta)}).simplify()
println_indented(sin2_eq)
theta_eq = Eq(theta, atan2(sin2_eq.rhs, cos2_eq.rhs) / 2)
println_indented(theta_eq)
println_indented(theta_eq.subs(ecc**2 * L, sign(ecc**2 * L)))
println_indented(theta_eq.subs(ecc**2 * L, 1 / sign(det)))
