#!/usr/bin/env python

from sympy import I, cancel, gcd, symbols

from lib.conic import conic_through_points, projective_conic_center
from lib.point import CIRCULAR_POINT_I, CIRCULAR_POINT_J
from research.poly_as_det import normalize_det3, poly_as_det3
from research.sympy_utils import println_indented

p1 = x1, y1, z1 = symbols("x1 y1 z1")
p2 = x2, y2, z2 = symbols("x2 y2 z2")
p3 = x3, y3, z3 = symbols("x3 y3 z3")

print("\nCircumcircle as conic matrix:\n")

circumcircle = conic_through_points(p1, p2, p3, CIRCULAR_POINT_I, CIRCULAR_POINT_J)
circumcircle /= -2 * I

println_indented(circumcircle.applyfunc(lambda e: normalize_det3(poly_as_det3(e))))

print("Circumcenter:\n")

circumcenter = projective_conic_center(circumcircle)
circumcenter = (circumcenter / gcd(list(circumcenter))).applyfunc(cancel)

println_indented(circumcenter.applyfunc(lambda e: normalize_det3(poly_as_det3(e))))
