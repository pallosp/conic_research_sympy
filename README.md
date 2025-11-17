# Calculating with conics (research project)

This project provides a symbolic framework for working with classical Euclidean
and projective geometry in 2D. Built on top of [SymPy](https://www.sympy.org/),
it offers tools to construct, transform, and analyze geometric entities such as
points, lines, and conics in projective space. It adapts several algorithms from
Jürgen Richter-Gebert’s Projective Geometry book. The framework is useful for
both symbolic and numeric computations, educational purposes, and geometric
reasoning.

The API supports:

- Conic representation in quadratic and polar form
- Constructing conics from five points, focus/directrix, center/radii, and more
- Line and point operations including bisectors, parallels, perpendiculars, etc.
- Geometric transformations like translation, scaling, and rotation in
  projective space
- Conic classification (e.g., ellipse, parabola, hyperbola, degenerate cases)
- Computation of conic properties (e.g. focus, eccentricity, asymptotes,
  vertices, ideal points)
- Incidence and distance calculations

[Complete API reference](docs/api.md)

## Development

### Installation

```sh
git clone https://github.com/pallosp/conic_research_sympy
poetry env use 3.13
poetry install
poetry run pre-commit install

# Optionally
`poetry env activate`
```

### Useful commands

Regenerating the API docs:

```sh
pydoc-markdown
```

Running the tests (with optional coverage report):

```sh
pytest

pytest --cov src/lib --cov-report html
open htmlcov/index.html
```

Running the linter (with optional auto-fixes):

```sh
ruff check
ruff check --fix
```

Running all pre-commit checks for the staged files:

```sh
pre-commit
```
