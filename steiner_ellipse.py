#!/usr/bin/env python

# Algorithm to construct an ellipse matrix from its center and 3 perimeter
# points.
#
# A special is the Steiner ellipse, when the ellipse's center coincides
# with the centroid of the triangle around the perimeter points.

from sympy import symbols, Matrix

from lib.central_conic import ConicCenter
from lib.matrix import QuadraticForm
from lib.transform import TransformConic, Translate

# Simplification WLOG: let the center point be at (0, 0).
#
# Equation of such ellipses:
#   a⋅x² + b⋅x⋅y + c⋅y² = 1
#
# The ellipse goes through (xᵢ,yᵢ) =>
#   Ɐi a⋅xᵢ² + b⋅xᵢ⋅yᵢ + c⋅yᵢ² = 1
#
# This linear equation system in matrix form:
#
#   [x₁⋅x₁ x₁⋅y₁ y₁⋅y₁] [a]   [1]
#   [x₂⋅x₂ x₂⋅y₂ y₂⋅y₂] [b] = [1]
#   [x₃⋅x₃ x₃⋅y₃ y₃⋅y₃] [c]   [1]
#
# Solution using Cramer's rule:
#
#   Mᵢ = M after replacing its ith column with [1, 1, 1]ᵀ
#   a = det M₁ / det M
#   b = det M₂ / det M
#   c = det M₃ / det M
#
# Ellipse formula:
#
#   det M₁⋅x² + det M₂⋅x⋅y + det M₃⋅y² - det M = 0
#
# In matrix form:
#
#   [det M₁    det M₂/2  0     ]
#   [det M₂/2  det M₃    0     ]
#   [0         0         -det M]

# Sanity check:

x1, y1, x2, y2, x3, y3 = symbols("x1,y1,x2,y2,x3,y3")

m = Matrix(
    [
        [x1 * x1, x1 * y1, y1 * y1],
        [x2 * x2, x2 * y2, y2 * y2],
        [x3 * x3, x3 * y3, y3 * y3],
    ]
)

m1, m2, m3 = m.copy(), m.copy(), m.copy()
m1[:, 0] = Matrix([1, 1, 1])
m2[:, 1] = Matrix([1, 1, 1])
m3[:, 2] = Matrix([1, 1, 1])

ellipse = Matrix(
    [
        [m1.det(), m2.det() / 2, 0],
        [m2.det() / 2, m3.det(), 0],
        [0, 0, -m.det()],
    ]
)

assert ConicCenter(ellipse) == (0, 0)
assert QuadraticForm(ellipse, Matrix([x1, y1, 1])).equals(0)
assert QuadraticForm(ellipse, Matrix([x2, y2, 1])).equals(0)
assert QuadraticForm(ellipse, Matrix([x3, y3, 1])).equals(0)

# Generalization for arbitrary center point (xc, yc)
#
#   1. Substitute xᵢ with xᵢ-xc and yᵢ with yᵢ-yc in the ellipse matrix.
#   2. Translate the result back with (xc, yc).
#
# The translation modifies the matrix elements as follows:
#
#   simplify(TransformConic(ConicMatrix(a, b, c, 0, 0, f), Translate(xc, yc))
#
#   [a b 0]    [a           b           -a⋅xc-b⋅yc           ]
#   [b c 0] -> [b           c           -b⋅xc-c⋅yc           ]
#   [0 0 f]    [-a⋅xc-b⋅yc  -b⋅xc-c⋅yc  f-a⋅xc²+b⋅xc⋅yc+c⋅yc²]
#
# The Steiner ellipse is a special case when xc = (x₁+x₂+x₃)/3 and
# yc = (y₁+y₂+y₃)/3.
#
# I couldn't find a way to express the last row and column of the matrix with a
# reasonably simple formula.

# Sanity check:

xc, yc = symbols("xc,yc")

substitutions = {
    x1: x1 - xc,
    y1: y1 - yc,
    x2: x2 - xc,
    y2: y2 - yc,
    x3: x3 - xc,
    y3: y3 - yc,
}

ellipse = TransformConic(ellipse.subs(substitutions), Translate(xc, yc))

assert ConicCenter(ellipse)[0].equals(xc)
assert QuadraticForm(ellipse, Matrix([x1, y1, 1])).equals(0)
