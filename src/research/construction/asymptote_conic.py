#!/usr/bin/env python

from typing import TYPE_CHECKING, cast

from sympy import gcd, pprint, simplify, symbols

from lib.central_conic import conic_center
from lib.conic import ideal_points
from lib.degenerate_conic import line_pair_conic
from lib.line import line_between
from lib.matrix import conic_matrix

if TYPE_CHECKING:
    from lib.type_aliases import PointVector

conic = conic_matrix(*symbols("a,b,c,d,e,f", positive=True))
# ideal_points() returns explicit vectors for this conic, which pyright can't
# know.
ideal_pt_1, ideal_pt_2 = cast("tuple[PointVector, PointVector]", ideal_points(conic))
center = conic_center(conic)
asymptote1 = line_between(center, ideal_pt_1)
asymptote2 = line_between(center, ideal_pt_2)
asymptote_conic = line_pair_conic(asymptote1, asymptote2)

asymptote_conic = asymptote_conic.applyfunc(simplify)
asymptote_conic /= gcd(asymptote_conic.flat())

print("\nDegenerate conic from the two asymptotes:\n")
pprint(asymptote_conic)
