#!/usr/bin/env python3

from collections.abc import Callable

from sympy import Expr, Matrix, gcd, pi, pprint, symbols

from lib.polar_conic import PointAtAngle

a, b, c, d, e, f, g, h, i, alpha, beta = symbols("a,b,c,d,e,f,g,h,i,alpha,beta")
polar_matrix = Matrix([[a, b, c], [d, e, f], [g, h, i]])


def GetIntersectionOfSecants(end_point_func: Callable[[Expr], Expr]) -> Matrix:
    alpha, beta = symbols("alpha,beta")

    p1 = PointAtAngle(polar_matrix, alpha)
    p2 = PointAtAngle(polar_matrix, end_point_func(alpha))
    p3 = PointAtAngle(polar_matrix, beta)
    p4 = PointAtAngle(polar_matrix, end_point_func(beta))

    secant1 = p1.cross(p2)
    secant2 = p3.cross(p4)

    common_point = secant1.cross(secant2)
    common_point /= gcd(common_point[0], common_point[1])
    for i in range(3):
        common_point[i] = common_point[i].factor()

    # For certain start point -> end point mapping functions (see below)
    # alpha and beta vanish.
    return common_point


print("\nPolar conic matrix:\n")
pprint(polar_matrix)

print("\nIntersections of secants between θ and π-θ:\n")
pprint(GetIntersectionOfSecants(lambda angle: pi - angle))

print("\nIntersections of secants between θ and -θ:\n")
pprint(GetIntersectionOfSecants(lambda angle: -angle) * -1)

print("\nIntersections of secants between θ and θ+π:\n")
pprint(GetIntersectionOfSecants(lambda angle: angle + pi))
