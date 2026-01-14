#!/usr/bin/env python

from sympy import gcd, pprint, simplify, symbols

from lib.central_conic import conic_center
from lib.conic import IdealPoints
from lib.degenerate_conic import line_pair_conic
from lib.line import line_between
from lib.matrix import conic_matrix

conic = conic_matrix(*symbols("a,b,c,d,e,f", positive=True))
ideal_point_1, ideal_point_2 = IdealPoints(conic)
center = conic_center(conic)
asymptote1 = line_between(center, ideal_point_1)
asymptote2 = line_between(center, ideal_point_2)
asymptote_conic = line_pair_conic(asymptote1, asymptote2)

asymptote_conic = asymptote_conic.applyfunc(simplify)
asymptote_conic /= gcd(list(asymptote_conic))

print("\nDegenerate conic from the two asymptotes:\n")
pprint(asymptote_conic)
