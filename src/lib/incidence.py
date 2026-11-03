from collections.abc import Callable, Sequence
from itertools import product

from sympy import Expr, Matrix, MatrixBase, MatrixExpr, expand
from sympy.core.logic import fuzzy_not

from lib.matrix import is_full_rank, quadratic_form, skew_matrix
from lib.point import point_to_vec3
from lib.type_aliases import Conic, Line, Point, PolarConic


def line_contains_point(
    line: Line | MatrixExpr,
    point: Point | MatrixExpr,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether `point` is on `line`.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if the result is
    undecidable. Accepts unevaluated `MatrixExpr` arguments.
    """
    # SymPy's stubs type MatrixExpr products as Expr (fix upstream).
    incidence = (line.T * point_to_vec3(point))[0, 0]  # pyright: ignore[reportIndexIssue]
    return simplifier(incidence).is_zero


def conic_contains_point(
    conic: Conic | MatrixExpr,
    point: Point | MatrixExpr,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Checks if a point lies on a conic.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if the result is
    undecidable. Accepts unevaluated `MatrixExpr` arguments.
    """
    return simplifier(quadratic_form(conic, point_to_vec3(point))).is_zero


def conic_contains_line(
    conic: Conic,
    line: Line,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Checks if a line lies on a conic.

    Takes an optional `simplifier` callback that simplifies the elements of the
    containment matrix before it gets compared to zero. Returns `None` if the
    result is undecidable.

    *Formula*:
    [research/incidence/conic_line_containment.py](../src/research/incidence/conic_line_containment.py)
    """
    skew = skew_matrix(line)
    return (skew * conic * skew).applyfunc(simplifier).is_zero_matrix


def is_tangent(
    conic: Conic,
    line: Line | MatrixExpr,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether a line is tangent to a conic, i.e. it intersects the conic
    at two coinciding points.

    Special cases:
     - Line pair conics: lines through the intersection of the pair are
       tangent, including the lines of the pair.
     - Double line conics: every line is tangent.

    Takes an optional `simplifier` callback that simplifies the tangency
    polynomial before it gets compared to zero. Returns `None` if the result is
    undecidable. Accepts unevaluated `MatrixExpr` lines.

    *Formula*: `lᵀ·adj(C)·l = 0`, i.e. the line is on the dual conic.
    """
    return simplifier(quadratic_form(conic.adjugate(), line)).is_zero


def is_inside(
    conic: Conic,
    point: Point,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether a point is inside a conic, i.e. no real tangent of the
    conic passes through it.

    Inside points are
     - the points within an ellipse, a parabola or a hyperbola branch;
     - all points for imaginary ellipses.

    Points on the conic are neither inside nor
    [outside](#incidence.is_outside). Neither are any points for
    [degenerate conics](#conic_classes.is_degenerate).

    Assumes a real conic and a real point. Takes an optional `simplifier`
    callback that simplifies the sign polynomial before its sign is determined.
    Returns `None` if undecidable.

    *Formula*: `det(C)·pᵀ·C·p > 0`<br>
    *Proof*:
    [research/incidence/point_inside_conic.py](../src/research/incidence/point_inside_conic.py)
    """
    return simplifier(_inside_sign_poly(conic, point)).is_positive


def is_outside(
    conic: Conic,
    point: Point,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether a point is outside a conic, i.e. two real tangents of the
    conic pass through it.

    Outside points are
     - the points outside an ellipse or a parabola;
     - the points between the branches of a hyperbola, including its center.

    Points on the conic are neither [inside](#incidence.is_inside) nor
    outside. Neither are any points for
    [degenerate conics](#conic_classes.is_degenerate).

    Assumes a real conic and a real point. Takes an optional `simplifier`
    callback that simplifies the sign polynomial before its sign is determined.
    Returns `None` if undecidable.

    *Formula*: `det(C)·pᵀ·C·p < 0`<br>
    *Proof*:
    [research/incidence/point_inside_conic.py](../src/research/incidence/point_inside_conic.py)
    """
    return simplifier(_inside_sign_poly(conic, point)).is_negative


def _inside_sign_poly(conic: Conic, point: Point) -> Expr:
    return conic.det() * quadratic_form(conic, point_to_vec3(point))


def polar_conic_contains_point(
    polar_conic: PolarConic,
    point: Point,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Checks if a point lies on a conic that is specified in polar form.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if the result is
    undecidable.

    *Formula*:
    [research/incidence/polar_conic_point_incidence.py](../src/research/incidence/polar_conic_point_incidence.py)
    """
    c, s, i = (polar_conic.adjugate() * point_to_vec3(point)).flat()
    return simplifier(c * c + s * s - i * i).is_zero


def are_collinear(
    points: Sequence[Point],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether n points are collinear.

    Takes an optional `simplifier` callback that simplifies the collinearity
    polynomial before it gets compared to zero. Returns `None` if undecidable.
    """
    if len(points) <= 2:
        return True
    # Convert the points to homogenous coordinates, and treat them as lines.
    lines = [point_to_vec3(p) for p in points]
    return are_concurrent(lines, simplifier=simplifier)


def are_concurrent(
    lines: Sequence[Line],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether n lines are concurrent, i.e. go through the same point.

    Takes an optional callback that simplifies the concurrence polynomial
    before it gets compared to zero. Returns `None` if undecidable.

    Leverages the projective point-line duality, and uses the collinearity
    formula described at
    <https://en.wikipedia.org/wiki/Incidence_(geometry)#Collinearity>
    """
    if len(lines) <= 2:
        return True
    lines_as_matrix = MatrixBase.hstack(*lines)
    if len(lines) == 3:
        return lines_as_matrix.det().is_zero
    return simplifier((lines_as_matrix * lines_as_matrix.T).det()).is_zero


def are_on_same_conic(
    points: Sequence[Point],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether up to 6 points lie on the same conic section.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if undecidable.

    *Formula*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 10.2 (Conics and Cross-Ratios)

    *Research*:
    [research/incidence/conconicity.py](../src/research/incidence/conconicity.py)
    """
    if len(points) < 6:
        return True
    if len(points) > 6:
        raise NotImplementedError
    p1, p2, p3, p4, p5, p6 = [point_to_vec3(p) for p in points]
    d = [
        MatrixBase.hstack(p, q, r).det()
        for p, q, r in product([p1, p2], [p3, p4], [p5, p6])
    ]
    incidence_poly = d[0] * d[3] * d[5] * d[6] - d[1] * d[2] * d[4] * d[7]
    return simplifier(incidence_poly).is_zero


def are_cocircular(
    points: Sequence[Point],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether n points lie on the same circle.

    Note that four collinear points, as well as three collinear points and an
    arbitrary ideal point are also considered cocircular.

    Takes an optional `simplifier` callback that simplifies the cocircularity
    determinant before it gets compared to zero. Returns `None` if undecidable.

    *Formula*: [Measuring cocircularity in a point set](
    https://upcommons.upc.edu/entities/publication/2001f976-4a50-4415-a706-79bbb2ca40bc)<br>
    *Own research*:
    [research/incidence/cocircularity.py](../src/research/incidence/cocircularity.py)
    """
    if len(points) < 4:
        return True

    cc_matrix_elements = []
    for p in points:
        x, y, z = point_to_vec3(p).flat()
        cc_matrix_elements.append([x * z, y * z, x**2 + y**2, z**2])

    cocircularity_matrix = Matrix(cc_matrix_elements)

    return fuzzy_not(is_full_rank(cocircularity_matrix, simplifier=simplifier))
