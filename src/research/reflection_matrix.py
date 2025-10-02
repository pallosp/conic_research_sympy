#!/usr/bin/env python3

from sympy import Expr, Matrix, factor, pprint, symbols

from lib.point import PerpendicularFoot
from lib.transform import TransformationFromSamples


def Reflect(point: tuple[Expr, Expr], axis: Matrix) -> tuple[Expr, Expr]:
    foot = PerpendicularFoot(point, axis)
    return (2 * foot[0] - point[0], 2 * foot[1] - point[1])


axis = Matrix(symbols("a b c"))
samples = [(1, 0), (0, 1), (-1, 0), (0, -1)]
transform = TransformationFromSamples(samples, [Reflect(s, axis) for s in samples])

print("\nTransformation matrix that reflects to the ax+by+c=0 line:\n")
pprint(transform.applyfunc(factor))
print()
