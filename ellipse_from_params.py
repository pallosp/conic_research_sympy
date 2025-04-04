#!/usr/bin/env python

from sympy import expand
from sympy import factor
from sympy import init_printing
from sympy import pprint
from sympy import simplify
from sympy import symbols

from lib.circle import UnitCircle
from lib.transform import Rotate
from lib.transform import TransformConic
from lib.transform import Translate
from lib.transform import ScaleXY

init_printing(use_unicode=True)


def main():
    cx, cy, r1, r2, angle = symbols("cx,cy,r1,r2,θ")
    unit_circle = UnitCircle()
    scaling = ScaleXY(r1, r2)
    rotation = Rotate(angle)
    translation = Translate(cx, cy)
    transformation = translation * rotation * scaling
    ellipse = TransformConic(unit_circle, transformation)

    a, b, c = symbols("a,b,c")
    for i in range(len(ellipse)):
        if i not in (0, 1, 3, 4):
            ellipse[i] = factor(ellipse[i], cx, cy)
            ellipse[i] = ellipse[i].subs(ellipse[0], "a")
            ellipse[i] = ellipse[i].subs(expand(ellipse[1]), "b")
            ellipse[i] = ellipse[i].subs(ellipse[4], "c")
            ellipse[i] = simplify(ellipse[i])
            ellipse[i] = factor(ellipse[i], cx, cy)
            ellipse[i] = ellipse[i].subs(expand(simplify(ellipse[1]) * 2), "b*2")

    print("Ellipse(cx, cy, r1, r2, θ):\n")
    pprint(ellipse)


main()
