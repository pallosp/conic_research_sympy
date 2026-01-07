#!/usr/bin/env python

from sympy import Eq, Matrix, simplify, symbols

from research.sympy_utils import println_indented

a, b, c = symbols("a,b,c", real=True)
min_eigen, max_eigen = symbols("min_eigen,max_eigen", real=True)

min_eigen_abc, max_eigen_abc = Matrix([[a, b], [b, c]]).eigenvals()

print("\nEigenvalues of [[a, b], [b, c]]:\n")

println_indented(Eq(min_eigen, min_eigen_abc))
println_indented(Eq(max_eigen, max_eigen_abc))
println_indented(Eq(min_eigen + max_eigen, min_eigen_abc + max_eigen_abc))
println_indented(Eq(max_eigen - min_eigen, max_eigen_abc - min_eigen_abc))
println_indented(Eq(max_eigen * min_eigen, simplify(max_eigen_abc * min_eigen_abc)))
