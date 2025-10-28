#!/usr/bin/env python

import itertools
import time
from collections.abc import Callable, Sequence

from sympy import Determinant, Expr, I, Matrix, S, factor, sqrt, symbols

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

    def save_poly(p: Expr) -> Expr:
        nonlocal poly
        poly = p
        return p

    are_on_same_conic([p1, p2, p3, p4, i, j], simplifier=save_poly)
    return poly


def cocircularity_matrix(p1: Matrix, p2: Matrix, p3: Matrix, p4: Matrix) -> Expr:
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


def measure_ms(func: Callable) -> float:
    start = time.time()
    func()
    return int((time.time() - start) * 10000) / 10


def eval_cocircularity_expr(
    expr: Expr,
    substitutions: Sequence[tuple[Expr, Expr]],
) -> float:
    _ = expr.subs(substitutions).doit().is_zero


points = [Matrix(symbols(f"x{i} y{i} z{i}")) for i in range(1, 5)]
coordinates = list(itertools.chain.from_iterable(points))
rational_substitutions = [(coordinates[i], S.One / i) for i in range(12)]
sqrt_substitutions = [(coordinates[i], 1 + sqrt(i)) for i in range(12)]


print("\nRichter-Gebert cocircularity polynomial:\n")

poly = richter_gebert_poly(*points)
print(f"  {poly}\n")

time_ms = measure_ms(lambda: richter_gebert_poly(*points))
print(f"  Creation: {time_ms} ms.")

time_ms = measure_ms(lambda: poly.subs(rational_substitutions).expand().is_zero)
print(f"  Evaluation (rational coordinates): {time_ms} ms.")

time_ms = measure_ms(lambda: poly.subs(sqrt_substitutions).expand().is_zero)
print(f"  Evaluation (radical coordinates): {time_ms} ms.")

time_ms = measure_ms(lambda: poly.expand().is_zero)
print(f"  Evaluation (symbolic coordinates): {time_ms} ms.")

print("\nCocircularity determinant:\n")

m = cocircularity_matrix(*points)
cocircularity_det = Determinant(m)
println_indented(cocircularity_det)
time_ms = measure_ms(lambda: cocircularity_matrix(*points))
print(f"  Creation: {time_ms} ms.")

time_ms = measure_ms(lambda: m.subs(rational_substitutions).det().is_zero)
print(f"  Evaluation (rational coordinates): {time_ms} ms.")

time_ms = measure_ms(lambda: m.subs(sqrt_substitutions).det().is_zero)
print(f"  Evaluation (radical coordinates): {time_ms} ms.")

time_ms = measure_ms(lambda: m.subs(sqrt_substitutions).det("laplace").expand().is_zero)
print(f"  Evaluation (radical coordinates, laplace): {time_ms} ms.")

time_ms = measure_ms(lambda: m.det("laplace").expand().is_zero)
print(f"  Evaluation (symbolic coordinates, laplace): {time_ms} ms.")

print("\nRatio of the two polynomials:\n")

println_indented(factor(poly.expand() / m.det()))
