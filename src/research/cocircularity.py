#!/usr/bin/env python

import time
from collections.abc import Sequence

from sympy import Determinant, Expr, I, Matrix, factor, symbols

from lib.incidence import are_on_same_conic
from lib.point import point_to_xy
from research.util import println_indented


def richter_gebert_poly(p1: Matrix, p2: Matrix, p3: Matrix, p4: Matrix) -> Expr:
    """Tells whether p1, p2, p3 and p4 are cocircular.

    Uses the formula derived at Jürgen Richter-Gebert, Perspectives on
    Projective Geometry, section 18.2 (Cocircularity).
    """
    i, j = Matrix([-I, 1, 0]), Matrix([I, 1, 0])
    poly = None

    def save_and_expand_poly(p: Expr) -> Expr:
        nonlocal poly
        poly = p
        return p.expand()

    are_on_same_conic([p1, p2, p3, p4, i, j], simplifier=save_and_expand_poly)
    return poly


def concircularity_matrix(p1: Matrix, p2: Matrix, p3: Matrix, p4: Matrix) -> Expr:
    """Tells whether p1, p2, p3 and p4 are cocircular.

    Uses the cocircularity determinant described in the paper

    XX Spanish Meeting on Computational Geometry, Santiago de Compostela
    July 3-5, 2023
    Measuring cocircularity in a point set

    https://egc23.web.uah.es/wp-content/uploads/2023/06/EGC23_paper_20.pdf
    """
    x1, y1 = point_to_xy(p1)
    x2, y2 = point_to_xy(p2)
    x3, y3 = point_to_xy(p3)
    x4, y4 = point_to_xy(p4)

    return Matrix(
        [
            [x1, y1, x1**2 + y1**2, 1],
            [x2, y2, x2**2 + y2**2, 1],
            [x3, y3, x3**2 + y3**2, 1],
            [x4, y4, x4**2 + y4**2, 1],
        ],
    )


def benchmark_cocircularity_expr_ms(expr: Expr, coordinates: Sequence[Expr]) -> int:
    start = time.time()
    expr = expr.subs(
        zip(
            coordinates,
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 1),
            strict=True,
        ),
    )
    _ = expr.doit().is_zero
    return int((time.time() - start) * 1000)


p1 = Matrix(symbols("x1 y1 z1"))
p2 = Matrix(symbols("x2 y2 z2"))
p3 = Matrix(symbols("x3 y3 z3"))
p4 = Matrix(symbols("x4 y4 z4"))
coordinates = list(p1) + list(p2) + list(p3) + list(p4)

print("\nRichter-Gebert cocircularity polynomial:\n")

poly = richter_gebert_poly(p1, p2, p3, p4)
println_indented(poly)
time_ms = benchmark_cocircularity_expr_ms(poly, coordinates)
print(f"  Its evaluation took {time_ms}ms.\n")

print("Cocircularity determinant:\n")

m = concircularity_matrix(p1, p2, p3, p4)
cocircularity_det = Determinant(m)
println_indented(cocircularity_det)
time_ms = benchmark_cocircularity_expr_ms(cocircularity_det, coordinates)
print(f"  Its evaluation took {time_ms}ms.\n")

print("Their ratio:\n")

println_indented(factor(poly.expand() / m.det()))
