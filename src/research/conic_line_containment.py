#!/usr/bin/env python

"""A line lies on a conic iff all of the following hold:
 - The conic is degenerate, i.e. it factors into two real or complex conjugate
   lines.
 - It's concurrent to the conic's lines, i.e. the pole point corresponding to
   the line is the zero vector.
 - At least two different points of the line lie on the conic.
   (0, -c, b), (c, 0, -a) and (-b, a, 0) are three points on the line, among
   which at least two are non-zero and different.

It's also possible to express this relationship with a single matrix equation.

Notations:
 - C: conic matrix
 - l: line vector
 - S: skew-symmetric matrix of the line vector

Skew-symmetric matrix identities:
 - Sᵀ = -S
 - Expresses cross products as matrix multiplication: l ⨯ v = S v
 - Its image, {S x | x ∈ ℝ³} or {l ⨯ x | x ∈ ℝ³} is l, because l ⨯ x · x = 0

The following statements are all equivalent:
 - l lines on the conic
 - ∀p ∈ l: pᵀ C p = 0
 - ∀x ∈ ℝ³: (S x)ᵀ C (S x) = 0
 - ∀x ∈ ℝ³: xᵀ (Sᵀ C S) x = 0
 - Sᵀ C S = 0
 - (-S) C S = 0
 - S C S = 0
"""

from sympy import Matrix, symbols

from lib.conic import pole_point
from lib.matrix import conic_matrix, quadratic_form, skew_matrix
from research.sympy_utils import print_indented, println_indented

conic = conic_matrix(*symbols("a,b,c,d,e,f"))
line = Matrix(symbols("x,y,z"))

print("\nConic matrix:\n")
println_indented(conic)

print("\nLine vector:\n")
println_indented(line)

skew = skew_matrix(line)
print("\nThree points on line are the columns of\n")
println_indented(skew)

print("\nCorresponding quadratic forms:")
print("(they must be zero if the line lies on the conic)\n")
for i in range(skew.cols):
    print_indented(quadratic_form(conic, skew.col(i)).expand())

print("\nskew * conic * skewᵀ:")
print("(notice that its diagonal elements are the quadratic forms)\n")
println_indented((skew * conic * skew.T).expand())

print("\nPole point corresponding to the line:")
print("(must be the zero vector if the line is concurrent to the conic's lines)\n")
println_indented(pole_point(conic, line).expand())
