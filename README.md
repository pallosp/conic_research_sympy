# Calculating with conics (research project)

This project provides a symbolic framework for working with classic Euclidean
and projective geometry in 2D. Built on top of [SymPy](https://www.sympy.org/),
it offers tools to construct, transform, and analyze geometric entities such as
points, lines, and conics in projective space.

The API supports:

- Conic representation in quadratic and polar forms
- Constructing conics from 5 points, focus/directrix, center/radii, and more
- Line and point operations including bisectors, parallels, perpendiculars, etc.
- Geometric transformations like translation, scaling, and rotation in
  projective space
- Conic classification (e.g., ellipse, parabola, hyperbola, degenerate)
- Distance computations

The library is based on foundational projective geometry principles, with
several algorithms adapted from Jürgen Richter-Gebert’s Projective Geometry
book. It's useful for symbolic computations, educational purposes, and geometric
reasoning.

[Complete API reference](docs/api.md)

## Development

### Installation

```sh
git clone https://github.com/pallosp/conic_research_sympy
python3 -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
```

### Documentation

To regenerate the API docs, run

```sh
pip install -e '.[docs]'
pydoc-markdown
```
