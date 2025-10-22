#!/usr/bin/env python

from sympy import Abs, Eq, Expr, Matrix, symbols

from lib.conic import conic_from_poly
from lib.distance import point_point_distance
from lib.sympy_utils import mul_eq, sub_eq
from research.util import println_indented

cx, cy = symbols("cx cy", real=True)  # center
focus = symbols("fx fy")
dx, dy = symbols("dx dy", real=True)  # center -> focus vector
focus1 = (cx + dx, cy + dy)
focus2 = (cx - dx, cy - dy)
radius = symbols("r", real=True)  # center-vertex distance


def ComputeConicFromFociAndRadius(
    focus1: tuple[Expr, Expr],
    focus2: tuple[Expr, Expr],
    radius: Expr,
    conic_type: str,
) -> Matrix:
    assert conic_type in ("Ellipse", "Hyperbola")

    print(f"{conic_type} equation from focal points and radius:\n")

    d1, d2 = symbols("d1 d2", real=True)
    conic_eq = Eq(
        radius,
        (d1 + d2) / 2 if conic_type == "Ellipse" else Abs(d1 - d2) / 2,
    )
    println_indented(conic_eq)

    conic_eq = mul_eq(conic_eq, conic_eq)
    println_indented(conic_eq)

    conic_eq = sub_eq(conic_eq, (d1**2 + d2**2) / 4)
    conic_eq = Eq(conic_eq.lhs, conic_eq.rhs.simplify())
    println_indented(conic_eq)

    conic_eq = mul_eq(conic_eq, conic_eq)
    conic_eq = sub_eq(conic_eq, conic_eq.rhs)
    println_indented(conic_eq)

    print("Substitute d₁ and d₂ with the actual distances, and collect x and y.\n")

    x, y = symbols("x y")
    point_on_hyperbola = (x, y)
    dist_focus1 = point_point_distance(focus1, point_on_hyperbola)
    dist_focus2 = point_point_distance(focus2, point_on_hyperbola)
    conic_eq = conic_eq.subs(d1, dist_focus1).subs(d2, dist_focus2)
    conic_eq = Eq(conic_eq.lhs.expand().collect((x, y)), 0)
    println_indented(conic_eq)

    print("In matrix form:\n")

    conic_matrix = conic_from_poly(conic_eq.lhs)
    println_indented(conic_matrix)
    return conic_matrix


ellipse = ComputeConicFromFociAndRadius(focus1, focus2, radius, "Ellipse")

print(f"{'-' * 135}\n")

hyperbola = ComputeConicFromFociAndRadius(focus1, focus2, radius, "Hyperbola")

print(f"{'-' * 135}\n")

print("They are the same:\n")

println_indented(
    Eq(
        symbols("ellipse-matrix") - symbols("hyperbola-matrix"),
        ellipse - hyperbola,
        evaluate=False,
    ),
)

assert ellipse == hyperbola
