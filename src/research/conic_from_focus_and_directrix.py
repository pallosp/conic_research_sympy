#!/usr/bin/env python

from sympy import Eq, MatAdd, MatMul, Matrix, gcd, symbols

from lib.distance import PointLineDistance, PointPointDistance
from research.util import println_indented

print("\nConic equation from focus, directrix and eccentricity:\n")

# eccentricity = |distance(p, focus) / distance(p, directrix)|

directrix = Matrix(symbols("a,b,c"))
a, b, _ = directrix
focus = symbols("fx,fy")
eccentricity = symbols("e")
point_on_conic = symbols("x,y")

distance_from_focus = PointPointDistance(point_on_conic, focus)
distance_from_directrix = PointLineDistance(point_on_conic, directrix)

conic_eq = (eccentricity * distance_from_directrix) ** 2 - distance_from_focus**2
conic_eq = (conic_eq * (a * a + b * b)).simplify()

println_indented(Eq(conic_eq, 0))


print("\nExpanded to matrix:\n")

x, y = point_on_conic
coeffs = conic_eq.expand().collect([x, y, x * y], evaluate=False)
conic_matrix = Matrix(
    [
        [coeffs[x * x], coeffs[x * y] / 2, coeffs[x] / 2],
        [coeffs[x * y] / 2, coeffs[y * y], coeffs[y] / 2],
        [coeffs[x] / 2, coeffs[y] / 2, coeffs[1]],
    ],
)
println_indented(conic_matrix)


print("\nAs a sum of two matrices:\n")

ecc_square = eccentricity**2
ecc_matrix = conic_matrix.applyfunc(
    lambda el: el.collect(eccentricity).coeff(ecc_square),
)

remainder = conic_matrix - ecc_square * ecc_matrix
remainder_gcd = gcd(list(remainder))
remainder = MatMul(
    remainder.applyfunc(lambda el: (el / remainder_gcd).factor()), remainder_gcd,
)

conic_eq = MatAdd(MatMul(ecc_matrix, ecc_square), remainder)
println_indented(conic_eq)


print("\nDeterminant of the conic matrix above:\n")

println_indented(conic_matrix.det().factor())


print("\nThe determinant is positive if")
print()
print("  • eccentricity > 0 (not a circle)")
print("  • a²+b² > 0 (the directrix is a finite line)")
print("  • ax+by+c ≠ 0 (the focus is not on the directix)")


print("\nIts discriminant:\n")

println_indented(conic_matrix[:2, :2].det().factor())
