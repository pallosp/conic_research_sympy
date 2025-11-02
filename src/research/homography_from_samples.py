#!/usr/bin/env python

from collections.abc import Sequence
from itertools import chain, combinations

from sympy import Determinant, Expr, MatMul, Matrix, S, factor, gcd, pi, symbols

from lib.transform import homography_from_samples, rotate, scale
from research.util import print_indented

"""Calculates the transformation that maps (xᵢ,yᵢ) to (uᵢ,vᵢ) for i=0..3.

Source: https://franklinta.com/2014/09/08/computing-css-matrix3d-transforms/
"""

x0, y0, x1, y1, x2, y2, x3, y3 = symbols("x0,y0,x1,y1,x2,y2,x3,y3")
u0, v0, u1, v1, u2, v2, u3, v3 = symbols("u0,v0,u1,v1,u2,v2,u3,v3")

# According to
# https://franklinta.com/2014/09/08/computing-css-matrix3d-transforms/
#
# we can get the coefficients of the transformation matrix
#
#   [a b c]
#   [d e f]
#   [g h 1]
#
# by solving the
#
#   m * [a,b,c,d,e,f,g,h]ᵀ = [u0,v0,u1,v1,u2,v2,u3,v3]ᵀ
#
# linear equation system where m is a 8x8 matrix below.


def transformation_from_2d_samples(
    source_points: Sequence[tuple[Expr, Expr]],
    target_points: Sequence[tuple[Expr, Expr]],
) -> Matrix:
    (x0, y0), (x1, y1), (x2, y2), (x3, y3) = source_points
    (u0, v0), (u1, v1), (u2, v2), (u3, v3) = target_points

    m = Matrix(
        [
            [x0, y0, 1, 0, 0, 0, -u0 * x0, -u0 * y0],
            [0, 0, 0, x0, y0, 1, -v0 * x0, -v0 * y0],
            [x1, y1, 1, 0, 0, 0, -u1 * x1, -u1 * y1],
            [0, 0, 0, x1, y1, 1, -v1 * x1, -v1 * y1],
            [x2, y2, 1, 0, 0, 0, -u2 * x2, -u2 * y2],
            [0, 0, 0, x2, y2, 1, -v2 * x2, -v2 * y2],
            [x3, y3, 1, 0, 0, 0, -u3 * x3, -u3 * y3],
            [0, 0, 0, x3, y3, 1, -v3 * x3, -v3 * y3],
        ],
    )
    uv = Matrix([u0, v0, u1, v1, u2, v2, u3, v3])
    coefficients = [*list(m.inv() * uv), S.One]
    return Matrix(3, 3, coefficients).applyfunc(factor)


################################################################################

# The general solution is quite big. Solve it for the special case where
# (xᵢ,yᵢ) = (±1, ±1)

print("\nTransformation to map the (±1, ±1) square to the (uᵢ,vᵢ) quadrilateral.")

horizontal_square = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
target = [(u0, v0), (u1, v1), (u2, v2), (u3, v3)]
transform = transformation_from_2d_samples(horizontal_square, target)

# Sanity checks

scale2 = transform.subs({u0: 2, v0: 2, u1: -2, v1: 2, u2: -2, v2: -2, u3: 2, v3: -2})
assert scale2 == scale(2)

rot90 = transform.subs({u0: -1, v0: 1, u1: -1, v1: -1, u2: 1, v2: -1, u3: 1, v3: 1})
assert rot90 == rotate(pi / 2)

# The transformation matrix

print("Elements of the normalized transformation matrix:\n")

transform *= transform[0].as_numer_denom()[1]
for el in transform:
    print_indented(el)

################################################################################

print("\nTransformation to map the |x|+|y|=1 square to the (uᵢ,vᵢ) quadrilateral.")

diamond = [(1, 0), (0, 1), (-1, 0), (0, -1)]
transform = transformation_from_2d_samples(diamond, target)

print("Elements of the normalized transformation matrix:\n")

transform *= transform[0].as_numer_denom()[1]
for el in transform:
    print_indented(el)

################################################################################

print(
    "\nElements of the transformation matrix "
    "mapping the (±1, ±1) square to (uᵢ,vᵢ,wᵢ):\n",
)

w0, w1, w2, w3 = symbols("w0 w1 w2 w3")
homogeneous_target = [
    (u0 / w0, v0 / w0),
    (u1 / w1, v1 / w1),
    (u2 / w2, v2 / w2),
    (u3 / w3, v3 / w3),
]
t = transformation_from_2d_samples(horizontal_square, homogeneous_target)
t /= gcd(t[0], t[1])

collectibles = [
    a * b
    for a, b in chain(
        combinations([u0, u1, u2, u3], 2),
        combinations([v0, v1, v2, v3], 2),
        combinations([w0, w1, w2, w3], 2),
    )
]

for el in t:
    print_indented(el.collect(collectibles))

################################################################################

print("\nExpressed with determinants:\n")

points = Matrix([[u0, u1, u2, u3], [v0, v1, v2, v3], [w0, w1, w2, w3]])
dets = [Determinant(points[:, [j for j in range(4) if j != i]]) for i in range(4)]
det_values = [det.doit() for det in dets]

left = points[:, :3]
right = (left.inv() * t).applyfunc(factor)
right = right.applyfunc(lambda el: el.subs(zip(det_values, dets, strict=True)))

print_indented(MatMul(left, right))

################################################################################

print("\nTransformation from (1,0,0), (0,1,0), (0,0,1), (1,1,1) to (uᵢ,vᵢ,wᵢ):\n")

special_source_points = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1))
t = homography_from_samples(special_source_points, homogeneous_target)
t *= w0 * w1 * w2 * w3
for el in t:
    print_indented(el)
