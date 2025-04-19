from sympy import Matrix, symbols

from lib.conic import ConicFromFocusAndDirectrix, Eccentricity


def test_eccentricity():
    a, b, c, fx, fy, e = symbols("a,b,c,fx,fy,e", nonnegative=True)
    conic = ConicFromFocusAndDirectrix((fx, fy), Matrix([a, b, c]), e)
    assert e == Eccentricity(conic).simplify()
