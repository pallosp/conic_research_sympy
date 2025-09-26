#!/usr/bin/env python

from sympy import gcd, pprint, simplify, symbols

from lib.central_conic import ConicCenter
from lib.conic import IdealPoints
from lib.degenerate_conic import LinePair
from lib.line import LineBetween
from lib.matrix import ConicMatrix

conic = ConicMatrix(*symbols("a,b,c,d,e,f", positive=True))
ideal_point_1, ideal_point_2 = IdealPoints(conic)
center = ConicCenter(conic)
asymptote1 = LineBetween(center, ideal_point_1)
asymptote2 = LineBetween(center, ideal_point_2)
asymptote_conic = LinePair(asymptote1, asymptote2)

asymptote_conic = asymptote_conic.applyfunc(simplify)
asymptote_conic /= gcd(list(asymptote_conic))

print("\nDegenerate conic from the two asymptotes:\n")
pprint(asymptote_conic)
