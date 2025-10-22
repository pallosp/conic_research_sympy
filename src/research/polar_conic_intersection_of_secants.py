#!/usr/bin/env python3

from collections.abc import Callable

from sympy import Expr, Matrix, factor, gcd, pi, symbols

from lib.polar_conic import point_at_angle
from research.util import println_indented

a, b, c, d, e, f, g, h, i, alpha, beta = symbols("a,b,c,d,e,f,g,h,i,alpha,beta")
polar_matrix = Matrix([[a, b, c], [d, e, f], [g, h, i]])


def get_intersection_of_secants(end_point_func: Callable[[Expr], Expr]) -> Matrix:
    alpha, beta = symbols("alpha,beta")

    p1 = point_at_angle(polar_matrix, alpha)
    p2 = point_at_angle(polar_matrix, end_point_func(alpha))
    p3 = point_at_angle(polar_matrix, beta)
    p4 = point_at_angle(polar_matrix, end_point_func(beta))

    secant1 = p1.cross(p2)
    secant2 = p3.cross(p4)

    common_point = secant1.cross(secant2)
    common_point /= gcd(common_point[0], common_point[1])

    # For certain start point -> end point mapping functions (see below)
    # alpha and beta vanish.
    return common_point.applyfunc(factor)


print("\nPolar conic matrix:\n")
println_indented(polar_matrix)

print("Intersections of secants between θ and π-θ:\n")
println_indented(get_intersection_of_secants(lambda angle: pi - angle))

print("Intersections of secants between θ and -θ:\n")
println_indented(get_intersection_of_secants(lambda angle: -angle) * -1)

print("Intersections of secants between θ and θ+π:\n")
println_indented(get_intersection_of_secants(lambda angle: angle + pi))
