#!/usr/bin/env python

import time
from collections.abc import Callable, Sequence

from sympy import Expr, Matrix, S, sqrt, symbols


def measure_det(
    matrix_size: int,
    element_gen: Callable[[int], Expr],
    method: str | None,
    *,
    expand: bool = False,
) -> float:
    elements = [element_gen(i) for i in range(matrix_size**2)]
    m = Matrix(matrix_size, matrix_size, elements)
    start_sec = time.time()
    det = m.det(method=method) if method else m.det()
    if expand:
        det.expand()
    return time.time() - start_sec


def benchmark_group(
    matrix_sizes: Sequence[int],
    element_gen: Callable[[int], Expr],
    matrix_type: str,
    method: str | None,
    *,
    expand: bool = False,
) -> None:
    for n in matrix_sizes:
        elapsed_sec = measure_det(n, element_gen, method, expand=expand)
        elapsed_ms = "%.1f" % (elapsed_sec * 1000)
        method_str = (method or "default") + " method"
        if expand:
            method_str += ", expanded"
        print(f"{n}x{n} {matrix_type} matrix, {method_str}: {elapsed_ms} ms")


methods = [None, "bareiss", "berkowitz", "bird", "domain-ge", "laplace", "lu"]

for method in methods:
    sizes = [4, 6, 8, 10]
    benchmark_group(sizes, lambda i: 1234 // (i + 1), "int", method, expand=False)
    print()

for method in [m for m in methods if m != "domain-ge"]:
    benchmark_group([6], lambda i: S.One / (i + 1), "rational", method, expand=False)
print()

for method in [m for m in methods if m != "domain-ge"]:
    benchmark_group([4], lambda i: 1 + sqrt(i), "1+√i", method, expand=False)
print()

# The default bareiss algorithm is ruled out for general slowness,
# and lu is ruled out for expansion slowness.
for method in ["berkowitz", "bird", "domain-ge", "laplace"]:
    benchmark_group([5], lambda i: symbols(f"x{i}"), "xᵢ", method, expand=True)
    benchmark_group([5], lambda i: symbols(f"x{i}") + 1, "xᵢ+1", method, expand=False)
    benchmark_group([5], lambda i: symbols(f"x{i}") + 1, "xᵢ+1", method, expand=True)
    print()
