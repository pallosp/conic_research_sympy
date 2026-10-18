from typing import TypeAlias

from sympy import Expr, MatrixBase

#: A scalar value: a SymPy expression or a Python number. Functions convert
#: numbers to SymPy expressions as needed.
Scalar: TypeAlias = Expr | int | float

#: A point given by its Euclidean coordinates `(x, y)` or its homogeneous
#: coordinates `(x, y, z)`, as a tuple, a list or a column vector.
Point: TypeAlias = (
    MatrixBase
    | tuple[Scalar, Scalar]
    | tuple[Scalar, Scalar, Scalar]
    # Lists are invariant, so list[Scalar] alone would reject e.g. list[Expr].
    | list[Expr]
    | list[int]
    | list[float]
    | list[Scalar]
)

#: Same as [Point](#type_aliases.Point), but conveys that the point is
#: expected to be finite, i.e. not an ideal point.
FinitePoint: TypeAlias = Point

#: A point as a homogeneous 3D column vector `(x, y, z)ᵀ`. Functions that
#: return projective points use this type; those returning Euclidean
#: coordinates return a 2D `Matrix` instead.
PointVector: TypeAlias = MatrixBase

#: A projective line `ax + by + c = 0`, represented as the 3D column vector
#: `(a, b, c)ᵀ`.
Line: TypeAlias = MatrixBase
