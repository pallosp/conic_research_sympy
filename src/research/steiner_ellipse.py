#!/usr/bin/env python

# Algorithm to construct a conic matrix from its center and three perimeter
# points.
#
# A special case is the Steiner ellipse, where the conic's center coincides
# with the centroid of the triangle formed by the perimeter points.

from sympy import Determinant, Eq, Expr, Matrix, symbols

from lib.central_conic import ConicCenter
from lib.incidence import ConicContainsPoint
from lib.transform import TransformConic, Translate
from research.util import println_indented

# Simplification WLOG: let the center point be at (0, 0).
#
# Equation of such conics:
#   a⋅x² + b⋅x⋅y + c⋅y² = 1
#
# The conic goes through (xᵢ,yᵢ) =>
#   Ɐi a⋅xᵢ² + b⋅xᵢ⋅yᵢ + c⋅yᵢ² = 1
#
# This linear equation system in matrix form:
#
#   [x₁² x₁⋅y₁ y₁²] [a]   [1]
#   [x₂² x₂⋅y₂ y₂²] [b] = [1]
#   [x₃² x₃⋅y₃ y₃²] [c]   [1]
#
# Solution using Cramer's rule:
#
#   Mᵢ = M after replacing its ith column with [1, 1, 1]ᵀ
#   a = det M₁ / det M
#   b = det M₂ / det M
#   c = det M₃ / det M
#
# Conic formula:
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
    ],
)

m1, m2, m3 = m.copy(), m.copy(), m.copy()
m1[:, 0] = Matrix([1, 1, 1])
m2[:, 1] = Matrix([1, 1, 1])
m3[:, 2] = Matrix([1, 1, 1])

conic = Matrix(
    [
        [m1.det(), m2.det() / 2, 0],
        [m2.det() / 2, m3.det(), 0],
        [0, 0, -m.det()],
    ],
)

assert ConicCenter(conic) == (0, 0)
assert ConicContainsPoint(conic, (x1, y1))
assert ConicContainsPoint(conic, (x2, y2))
assert ConicContainsPoint(conic, (x3, y3))

# Generalization for arbitrary center point (xc, yc)
#
#   1. Substitute xᵢ with xᵢ-xc and yᵢ with yᵢ-yc in the conic matrix.
#   2. Translate the result back with (xc, yc).
#
# The translation modifies the matrix elements as follows:
#
#   simplify(TransformConic(ConicMatrix(a, b, c, 0, 0, f), Translate(xc, yc))
#
#   [a b 0]    [a             b             -a⋅xc - b⋅yc               ]
#   [b c 0] -> [b             c             -b⋅xc - c⋅yc               ]
#   [0 0 f]    [-a⋅xc - b⋅yc  -b⋅xc - c⋅yc  f - a⋅xc² + b⋅xc⋅yc + c⋅yc²]
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

conic = TransformConic(conic.subs(substitutions), Translate(xc, yc))

assert ConicCenter(conic)[0].equals(xc)
assert ConicContainsPoint(conic, Matrix([x1, y1, 1]))

# The Steiner ellipse is a special case when xc = (x₁+x₂+x₃)/3 and
# yc = (y₁+y₂+y₃)/3.

steiner = conic.subs({xc: (x1 + x2 + x3) / 3, yc: (y1 + y2 + y3) / 3})

# Simpify with the common factor in the matrix elements
steiner /= Matrix([[x1, x2, x3], [y1, y2, y3], [1, 1, 1]]).det() / -18

a, _, _, b, c, _, d, e, f = (elem.factor() for elem in steiner)


def PolyAsDet(
    poly: Expr,
    col1_options: list[Matrix],
    col2_options: list[Matrix],
    col3_options: list[Matrix],
) -> Determinant:
    """Expresses poly as the determinant of a 3x3 matrix."""
    for c1 in col1_options:
        for c2 in col2_options:
            for c3 in col3_options:
                matrix = Matrix.hstack(c1, c2, c3)
                if poly.equals(matrix.det()):
                    return Determinant(matrix)
    raise AssertionError("not found")


col_ones = Matrix([1, 1, 1])
col_x_options = [
    Matrix([x1, x2, x3]),
    Matrix([x2 - x3, x3 - x1, x1 - x2]),
    Matrix([x2 + x3, x3 + x1, x1 + x2]),
]
col_y_options = [
    Matrix([y1, y2, y3]),
    Matrix([y2 - y3, y3 - y1, y1 - y2]),
    Matrix([y2 + y3, y3 + y1, y1 + y2]),
]
col_xy_asymmetric = Matrix(
    [
        x1 * y2 - x2 * y1,
        x2 * y3 - x3 * y2,
        x3 * y1 - x1 * y3,
    ],
)
col_xy_symmetric = [
    Matrix(
        [
            x1 * (y2 - y3) - y1 * (x2 - x3),
            x2 * (y3 - y1) - y2 * (x3 - x1),
            x3 * (y1 - y2) - y3 * (x1 - x2),
        ],
    ),
    Matrix([x2 * y3 - x3 * y2, x3 * y1 - x1 * y3, x1 * y2 - x2 * y1]),
]


print("\nElements of the Steiner ellipse matrix:\n")

println_indented(
    Eq(symbols("a"), PolyAsDet(a, col_y_options, col_y_options, [col_ones])),
)
println_indented(
    Eq(symbols("b"), PolyAsDet(b, col_x_options, col_y_options, [col_ones])),
)
println_indented(
    Eq(symbols("c"), PolyAsDet(c, col_x_options, col_x_options, [col_ones])),
)
println_indented(
    Eq(symbols("d"), PolyAsDet(d, col_y_options, col_y_options, col_x_options)),
)
println_indented(
    Eq(symbols("e"), PolyAsDet(e, col_x_options, col_x_options, col_y_options)),
)
println_indented(
    Eq(symbols("f"), PolyAsDet(f, col_x_options, col_xy_symmetric, col_y_options)),
)
println_indented(
    Eq(
        symbols("f") / 2,
        PolyAsDet(f / 2, col_x_options, [col_xy_asymmetric], col_y_options),
    ),
)

# Steiner inellipse
#
# It's essentially the Steiner ellipse scaled to half size from its center.
# According to scale_conic_from_center.py its a, b, c, d and e coefficients are
# the same (modulo a constant factor).

print("\nElements of the Steiner inellipse matrix:\n")

inellipse_subs = {
    x1: (x2 + x3) / 2,
    x2: (x3 + x1) / 2,
    x3: (x1 + x2) / 2,
    y1: (y2 + y3) / 2,
    y2: (y3 + y1) / 2,
    y3: (y1 + y2) / 2,
}

a, b, c, d, e, f = (v.xreplace(inellipse_subs) * 4 for v in [a, b, c, d, e, f])

println_indented(
    Eq(symbols("a"), PolyAsDet(a, col_y_options, col_y_options, [col_ones])),
)
println_indented(
    Eq(symbols("b"), PolyAsDet(b, col_x_options, col_y_options, [col_ones])),
)
println_indented(
    Eq(symbols("c"), PolyAsDet(c, col_x_options, col_x_options, [col_ones])),
)
println_indented(
    Eq(symbols("d"), PolyAsDet(d, col_y_options, col_y_options, col_x_options)),
)
println_indented(
    Eq(symbols("e"), PolyAsDet(e, col_x_options, col_x_options, col_y_options)),
)
println_indented(
    Eq(
        symbols("f") * 2,
        PolyAsDet(f * 2, col_x_options, col_xy_symmetric, col_y_options),
    ),
)
