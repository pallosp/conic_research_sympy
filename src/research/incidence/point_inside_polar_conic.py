#!/usr/bin/env python

"""Decides whether a point is inside or outside a conic given in polar form.

Notations:
 - polar: real polar conic matrix, mapping the unit circle to the conic
 - unit_circle = diag(-1, -1, 1): the circle 1 - x² - y² = 0
 - conic: the same conic in quadratic form
 - point: real point vector
 - [u, v, w]ᵀ = adj(polar)·point: the point's preimage, up to scale

See research/incidence/point_inside_conic.py for the general criterion: the
point is inside the conic iff det(conic)·pointᵀ·conic·point > 0, and outside
iff it's < 0.
"""

from typing import cast

from sympy import Eq, Matrix, MatrixBase, Symbol, symbols

from lib.circle import UNIT_CIRCLE
from lib.matrix import quadratic_form
from lib.polar_conic import conic_from_polar_matrix
from research.sympy_utils import println_indented


def equals(lhs: MatrixBase, rhs: MatrixBase) -> bool | None:
    # SymPy's stubs lack MatrixBase.equals (fix upstream).
    return cast("Matrix", lhs).equals(rhs)


polar = Matrix(3, 3, symbols("a:i", real=True))
point = Matrix(symbols("x:z", real=True))
conic = conic_from_polar_matrix(polar)

print("\n1. The conic is the image of the unit circle under the polar matrix, so")
print("   conic = adj(polar)ᵀ·unit_circle·adj(polar).")
assert equals(conic, polar.adjugate().T * UNIT_CIRCLE * polar.adjugate())

print("\n2. Hence pointᵀ·conic·point")
print("     = (adj(polar)·point)ᵀ·unit_circle·(adj(polar)·point)")
print("     = w² - u² - v².")
u, v, w = (polar.adjugate() * point).flat()
assert (quadratic_form(conic, point) - (w * w - u * u - v * v)).expand() == 0

print("\n3. det(conic) = det(adj(polar))²·det(unit_circle) = det(polar)⁴.")
assert (conic.det() - polar.det() ** 4).expand() == 0

print("\n4. Therefore\n")
println_indented(
    Eq(
        Symbol("det(conic)·pointᵀ·conic·point"),
        Symbol("det(polar)⁴") * Symbol("(w² - u² - v²)"),
    ),
)

print("5. For an invertible polar matrix, det(polar)⁴ > 0, so the point is")
print("    - inside the conic iff w² > u² + v², i.e. its preimage is inside the")
print("      unit circle;")
print("    - outside iff w² < u² + v², i.e. its preimage is outside.")
print("   This agrees with projective transformations preserving tangency.")

print("\n6. For a singular polar matrix, det(conic) = 0, so the conic is")
print("   degenerate, and no point is inside or outside. w² - u² - v² may still")
print("   be non-zero, so the implementation checks det(polar) ≠ 0 separately.\n")
