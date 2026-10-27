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
#: coordinates return a 2D `Matrix` instead, except
#: [point_to_xy](#point.point_to_xy), which returns a tuple.
PointVector: TypeAlias = MatrixBase

#: A direction given by a 2D vector `(dx, dy)` or by the ideal point
#: `(dx, dy, 0)`, as a tuple, a list or a column vector. The vector's length
#: doesn't matter, only its direction.
Direction: TypeAlias = (
    MatrixBase
    | tuple[Scalar, Scalar]
    | tuple[Scalar, Scalar, Scalar]
    # Lists are invariant, so list[Scalar] alone would reject e.g. list[Expr].
    | list[Expr]
    | list[int]
    | list[float]
    | list[Scalar]
)

#: A projective line `ax + by + c = 0`, represented as the 3D column vector
#: `(a, b, c)ᵀ`.
Line: TypeAlias = MatrixBase

#: A conic section `ax² + 2bxy + cy² + 2dx + 2ey + f = 0`, represented as the
#: 3x3 symmetric matrix `[[a, b, d], [b, c, e], [d, e, f]]`.
Conic: TypeAlias = MatrixBase

#: A conic in polar parametric form, represented as a 3x3 matrix `M` such that
#: `M * (cos θ, sin θ, 1)ᵀ` traces the homogeneous coordinates of its points.
PolarConic: TypeAlias = MatrixBase

#: A projective transformation of the plane, represented as a 3x3 matrix that
#: maps the homogeneous coordinates of points.
Transformation: TypeAlias = MatrixBase
