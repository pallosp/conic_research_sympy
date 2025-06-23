#!/usr/bin/env python

from sympy import Matrix, pprint, symbols, sympify

"""Calculates the transformation that maps (xᵢ,yᵢ) to (uᵢ,vᵢ) for i=0..3.
ᵢ
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
    ]
)
uv = Matrix([u0, v0, u1, v1, u2, v2, u3, v3])

# The general solution is quite big. Solve it for the special case where
# (xᵢ,yᵢ) = (±1, ±1)

map = {x0: 1, y0: 1, x1: -1, y1: 1, x2: -1, y2: -1, x3: 1, y3: -1}
coefficients = list(m.subs(map).inv() * uv) + [sympify(1)]
for i in range(len(coefficients)):
    coefficients[i] = coefficients[i].factor()
transform = Matrix(3, 3, coefficients)

# Sanity checks

scale2 = transform.subs({u0: 2, v0: 2, u1: -2, v1: 2, u2: -2, v2: -2, u3: 2, v3: -2})
assert scale2 == Matrix.diag([2, 2, 1])

rot90 = transform.subs({u0: -1, v0: 1, u1: -1, v1: -1, u2: 1, v2: -1, u3: 1, v3: 1})
assert rot90 == Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

# The transformation matrix

transform *= transform[0].as_numer_denom()[1]
print("\nElements of the transformation matrix:\n")
for el in transform:
    pprint(el)

# Another interesting case is when the source points are (±1, 0) and (0, ±1),
# although the solution is less symmetric.
#
# map = {x0: 1, y0: 0, x1: 0, y1: 1, x2: -1, y2: 0, x3: 0, y3: -1}
