#!/usr/bin/env python

import itertools
import math
import time
from collections.abc import Sequence

from sympy import Determinant, Expr, Matrix, symbols

from lib.point import point_to_vec3
from research.sympy_utils import println_indented


def richter_gebert_poly(six_points: Matrix) -> Expr:
    """Tells whether six points lie on the same conic.

    Uses the formula derived at Jürgen Richter-Gebert, Perspectives on
    Projective Geometry, section 10.2 (Conics and Cross-Ratios)
    """
    d = []
    for p in six_points[0:2]:
        for q in six_points[2:4]:
            for r in six_points[4:6]:
                d.append(Matrix.hstack(p, q, r).det())  # noqa: PERF401
    return d[0] * d[3] * d[5] * d[6] - d[1] * d[2] * d[4] * d[7]


def conconicity_matrix(points: Matrix) -> Expr:
    """Tells whether six points lie on the same conic.

    Computes the conconicity determinant.
    """

    def to_matrix_row(point: Matrix) -> Matrix:
        x, y, z = point_to_vec3(point)
        return [x * x, y * y, z * z, x * y, y * z, z * x]

    return Matrix([to_matrix_row(p) for p in points])


def fibonacci(n: int) -> int:
    return round(((math.sqrt(5) + 1) / 2) ** n / math.sqrt(5))


def benchmark_conconicity_expr_ms(
    expr: Expr,
    substitutions: Sequence[tuple[Expr, Expr]],
) -> int:
    start = time.time()
    _ = expr.subs(substitutions).doit().is_zero
    return int((time.time() - start) * 1000)


points = [Matrix(symbols(f"x{i} y{i} z{i}")) for i in range(1, 7)]
coordinates = list(itertools.chain.from_iterable(points))
substitutions = [(coordinates[i], fibonacci(i)) for i in range(18)]

print("\nRichter-Gebert conconicity polynomial:\n")

poly = richter_gebert_poly(points)
print(f"  {poly}\n")
time_ms = benchmark_conconicity_expr_ms(poly, substitutions)
print(f"  Its evaluation for integer coordinates took {time_ms}ms.\n")

print("Conconicity determinant:\n")

m = conconicity_matrix(points)
conconicity_det = Determinant(m)
println_indented(conconicity_det)
time_ms = benchmark_conconicity_expr_ms(conconicity_det, substitutions)
print(f"  Its evaluation for integer coordinates took {time_ms}ms.\n")

print("Checking whether the two polynomials are the same:\n")

assert poly.subs(substitutions) == conconicity_det.subs(substitutions).doit()
print("  OK\n")
