#!/usr/bin/env python

from textwrap import indent
from sympy import Eq, Matrix, Piecewise, Pow, factor, pretty, sign, symbols

from lib.conic import ConicFromFocusAndDirectrix


def pprint_indented(expr):
    print(indent(pretty(expr), "  "))


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
    pprint_indented(eq)

print("\nNormalize the directrix so that a² + b² = 1:\n")

for i in range(len(coeff_eq)):
    coeff_eq[i] = coeff_eq[i].subs(a * a + b * b, 1)
    pprint_indented(coeff_eq[i])

print("\nDeterminant of the conic matrix:\n")

det = symbols("det")
det_eq = Eq(det, conic.det().factor().subs(a * a + b * b, 1))
pprint_indented(det_eq)

print("\nThis implies that sign(λ) = sign(det).")
print("Note that the formula breaks down for circles.\n")

print("Calculate the eigenvalues of the top left submatrix.\n")

p, q, r = symbols("p,q,r")
sorted_eigenvalues = [
    v.subs({p: conic[0], q: conic[3], r: conic[4]})
    .replace(Pow, lambda base, exp: Pow(base.factor(), exp))
    .factor()
    .subs(a * a + b * b, 1)
    .collect(ecc)
    for v in Matrix([[p, q], [q, r]]).eigenvals()
]


def Add(*eqs):
    left = eqs[0].lhs
    right = eqs[0].rhs
    for term in eqs[1:]:
        left += term.lhs
        right += term.rhs
    return Eq(left, right)


def Sub(eq0, eq1):
    return Eq(eq0.lhs - eq1.lhs, eq0.rhs - eq1.rhs)


def Div(eq, denom):
    return Eq(eq.lhs / denom, eq.rhs / denom)


def Swap(eq):
    return Eq(eq.rhs, eq.lhs)


min_eigen, max_eigen = symbols("min_eigen,max_eigen")
eigenvalue_min_eq = Eq(min_eigen, sorted_eigenvalues[0])
eigenvalue_max_eq = Eq(max_eigen, sorted_eigenvalues[1])
eigenvalue_sum_eq = factor(Add(eigenvalue_min_eq, eigenvalue_max_eq))
eigenvalue_diff_eq = factor(Sub(eigenvalue_max_eq, eigenvalue_min_eq))

eigenvalue_eqs = [
    eigenvalue_min_eq,
    eigenvalue_max_eq,
    eigenvalue_sum_eq,
    eigenvalue_diff_eq,
]

for eq in eigenvalue_eqs:
    pprint_indented(eq)

print("\nAs sign(λ) = sign(det):\n")

eigenvalue_diff_eq = Eq((max_eigen - min_eigen) * sign(det, evaluate=False), ecc**2 * L)
pprint_indented(eigenvalue_diff_eq)
print()
lambda_eq = factor(Div(Sub(eigenvalue_diff_eq, eigenvalue_sum_eq), 2))
pprint_indented(lambda_eq)
lambda_eq = Eq(L, Piecewise((-min_eigen, det > 0), (-max_eigen, True)))
pprint_indented(lambda_eq)
ecc_square_eq = Swap(Div(eigenvalue_diff_eq, L))
print()
pprint_indented(ecc_square_eq)
