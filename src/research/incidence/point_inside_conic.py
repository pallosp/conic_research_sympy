#!/usr/bin/env python

"""Decides whether a point is inside or outside a conic, i.e. whether real
tangents pass through it.

Notations:
 - C: real conic matrix
 - p: real point vector
 - T = skew(p)ᵀ·adj(C)·skew(p): the tangent pair conic through p
"""

from typing import cast

from sympy import Eq, I, Matrix, MatrixBase, Symbol, symbols

from lib.conic import tangent_pair_conic
from lib.matrix import conic_matrix, quadratic_form, skew_matrix
from research.sympy_utils import println_indented


def equals(lhs: MatrixBase, rhs: MatrixBase) -> bool | None:
    # SymPy's stubs lack MatrixBase.equals (fix upstream).
    return cast("Matrix", lhs).equals(rhs)


conic = conic_matrix(*symbols("a:f", real=True))
point = Matrix(symbols("x:z", real=True))

print("\n1. A line pair conic with intersection point q has adj = -q·qᵀ")
print("   (see research/matrices/line_pair_adjugate.py).")

print("\n2. q is real for real lines. For complex conjugate lines l = u + 𝑖·v")
print("   and l̄ = u - 𝑖·v, q is purely imaginary:\n")

u = Matrix(symbols("u:3", real=True))
v = Matrix(symbols("v:3", real=True))
q = (u + I * v).cross(u - I * v)
println_indented(q.expand())
assert equals(q, -2 * I * u.cross(v))

print("3. Hence T consists of real lines iff adj(T) = -q·qᵀ is negative")
print("   semidefinite, and of complex conjugate lines iff it's positive")
print("   semidefinite.")

print("\n4. adj(T) = adj(skew(p)) · adj(adj(C)) · adj(skew(p)ᵀ), where")
print("   adj(skew(p)) = adj(skew(p)ᵀ) = p·pᵀ:\n")

skew = skew_matrix(point)
println_indented(skew.adjugate())
assert equals(skew.adjugate(), point * point.T)
assert equals(skew.T.adjugate(), point * point.T)

print("   and adj(adj(C)) = det(C)·C")
assert equals(conic.adjugate().adjugate(), conic.det() * conic)

print("\n5. Therefore adj(T) = p·pᵀ · det(C)·C · p·pᵀ = det(C)·(pᵀ·C·p)·p·pᵀ\n")

adj = tangent_pair_conic(conic, point).adjugate()
sign_poly = conic.det() * quadratic_form(conic, point)
assert equals(adj, sign_poly * point * point.T)
println_indented(Eq(Symbol("det(C)·pᵀ·C·p"), sign_poly.factor()))

print("6. Since p·pᵀ is positive semidefinite, for a non-degenerate real conic")
print("   and a real point")
print("    - det(C)·pᵀ·C·p < 0: two real tangents, the point is outside;")
print("    - det(C)·pᵀ·C·p > 0: two complex conjugate tangents, it's inside;")
print("    - pᵀ·C·p = 0: the point is on the conic, and T is its tangent, twice.")

print("\n7. For degenerate conics det(C) = 0, so adj(T) = 0, i.e. T has rank ≤ 1.")
print("   Through any point there's either a single real tangent, the line to")
print("   the singular point, or every line is tangent. So no point is inside")
print("   or outside.\n")
