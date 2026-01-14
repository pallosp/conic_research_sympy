#!/usr/bin/env python

import time
from collections.abc import Callable, Sequence

from sympy import Expr, Matrix, S, sqrt, symbols

WARMUP_TIMES = 1


def measure_det(
    matrix_size: int,
    element_gen: Callable[[int], Expr],
    method: str | None,
    *,
    expand: bool = False,
) -> float:
    elements = [element_gen(i) for i in range(matrix_size**2)]
    m = Matrix(matrix_size, matrix_size, elements)
    for _ in range(WARMUP_TIMES):
        _ = m.det(method=method) if method else m.det()

    start_sec = time.time()
    det = m.det(method=method) if method else m.det()
    if expand:
        det.expand()
    return time.time() - start_sec


def measure_rank(matrix_size: int, element_gen: Callable[[int], Expr]) -> float:
    elements = [element_gen(i) for i in range(matrix_size**2)]
    m = Matrix(matrix_size, matrix_size, elements)
    for _ in range(WARMUP_TIMES):
        _ = m.rank()
    start_sec = time.time()
    _ = m.rank()
    return time.time() - start_sec


def benchmark_group(
    matrix_sizes: Sequence[int],
    element_gen: Callable[[int], Expr],
    matrix_type: str,
    subject: str,  # det | expanded_det | rank
    det_method: str | None = None,
) -> None:
    for n in matrix_sizes:
        elapsed_sec = None
        method_str = None
        if subject == "rank":
            elapsed_sec = measure_rank(n, element_gen)
            method_str = "rank"
        else:
            expand = subject == "expanded_det"
            elapsed_sec = measure_det(n, element_gen, det_method, expand=expand)
            method_str = (det_method or "default") + " method"
            if expand:
                method_str += ", expanded"

        elapsed_ms = "%.1f" % (elapsed_sec * 1000)

        print(f"{n}x{n} {matrix_type} matrix, {method_str}: {elapsed_ms} ms")


methods = [None, "bareiss", "berkowitz", "bird", "domain-ge", "laplace", "lu"]

for method in methods:
    sizes = [4, 6, 8, 10]
    benchmark_group(sizes, lambda i: 1234 // (i + 1), "int", "det", method)
    print()

sizes = [4, 6, 8, 10]
benchmark_group(sizes, lambda i: 1234 // (i + 1), "int", "rank")
print()

for method in [m for m in methods if m != "domain-ge"]:
    benchmark_group([6], lambda i: S.One / (i + 1), "rational", "det", method)
benchmark_group([6], lambda i: S.One / (i + 1), "rational", "rank")
print()

for method in [m for m in methods if m != "domain-ge"]:
    benchmark_group([4], lambda i: 1 + sqrt(i), "1+√i", "expanded_det", method)
print()

# The default bareiss algorithm is ruled out for general slowness,
# and lu is ruled out for expansion slowness.
for method in ["berkowitz", "bird", "domain-ge", "laplace"]:
    benchmark_group([5], lambda i: symbols(f"x{i}"), "xᵢ", "expanded_det", method)
    benchmark_group([5], lambda i: symbols(f"x{i}") + 1, "xᵢ+1", "det", method)
    benchmark_group([5], lambda i: symbols(f"x{i}") + 1, "xᵢ+1", "expanded_det", method)
    print()
