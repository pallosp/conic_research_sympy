#!/usr/bin/env python

from sympy import Matrix, Pow, pprint, simplify, symbols

a, b, c = symbols("a,b,c", real=True)

parabola = Matrix([[a, b], [b, c]])
axis_dirs = [simplify(r[2][0]) for r in parabola.eigenvects()]

# Simplify using the fact that a*c - b*b = 0
axis_dirs = [dir.subs(b**2, a * c) for dir in axis_dirs]

# For simplify the sqrt subexpressions
for i, dir in enumerate(axis_dirs):
    for coord in dir:
        for atom in coord.atoms(Pow):
            if atom.exp.equals(1 / 2):
                axis_dirs[i] = axis_dirs[i].subs(atom.base, atom.base.factor())

print("\nThe axis direction of the parabola is one of\n")
pprint(axis_dirs)
