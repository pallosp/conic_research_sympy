#!/usr/bin/env python

from sympy import Matrix, pprint, symbols

from lib.conic import PolePoint
from lib.matrix import ConicMatrix, QuadraticForm, SkewMatrix

conic = ConicMatrix(*symbols("a,b,c,d,e,f"))
line = Matrix(symbols("x,y,z"))

print("\nConic matrix:\n")
pprint(conic)

print("\nLine vector:\n")
pprint(line)

skew = SkewMatrix(line)
print("\nThree points on line are the columns of\n")
pprint(skew)

print("\nCorresponding quadratic forms:")
print("(they must be zero if the line lies on the conic)\n")
for i in range(skew.cols):
    pprint(QuadraticForm(conic, skew.col(i)).expand())

print("\nskew * conic * skewᵀ:")
print("(notice that its diagonal elements are the quadratic forms)\n")
pprint((skew * conic * skew.T).expand())

print("\nPole point corresponding to the line:")
print("(must be the zero vector if the line is concurrent to the conic's lines)\n")
pprint(PolePoint(conic, line).expand())

print(
    "\nTODO: Why is the pole point zeroness check equivalent to that of "
    "the non-diagonal elements of skew * conic * skewᵀ?\n",
)
