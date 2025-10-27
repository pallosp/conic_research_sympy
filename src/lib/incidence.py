from collections.abc import Callable, Sequence

from sympy import Expr, Matrix, expand

from lib.matrix import quadratic_form, skew_matrix
from lib.point import point_to_vec3


def line_contains_point(
    line: Matrix,
    point: Matrix | Sequence[Expr],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether `point` is on `line`.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if the result is
    undecidable.
    """
    return simplifier(line.dot(point_to_vec3(point))).is_zero


def conic_contains_point(
    conic: Matrix,
    point: Matrix | Sequence[Expr],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Checks if a point lies on a conic.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if the result is
    undecidable.
    """
    return simplifier(quadratic_form(conic, point_to_vec3(point))).is_zero


def conic_contains_line(
    conic: Matrix,
    line: Matrix,
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Checks if a line lies on a conic.

    Takes an optional `simplifier` callback that simplifies the elements of the
    containment matrix before it gets compared to zero. Returns `None` if the
    result is undecidable.

    *Formula*:
    [research/conic_line_containment.py](../src/research/conic_line_containment.py)
    """
    skew = skew_matrix(line)
    return (skew * conic * skew).applyfunc(simplifier).is_zero_matrix


def are_collinear(
    points: Sequence[Matrix],
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
    lines: Sequence[Matrix],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether n lines are concurrent, i.e. go through the same point.

    Takes an optional callback that simplifies the concurrence polynomial
    before it gets compared to zero. Returns `None` if undecidable.

    Leverages the projective point-line duality, and uses the collinearity
    formula described at
    https://en.wikipedia.org/wiki/Incidence_(geometry)#Collinearity
    """
    if len(lines) <= 2:
        return True
    lines_as_matrix = Matrix.hstack(*lines)
    if len(lines) == 3:
        return lines_as_matrix.det().is_zero
    return simplifier((lines_as_matrix * lines_as_matrix.T).det()).is_zero


def are_on_same_conic(
    points: Sequence[Matrix],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether up to 6 points lie on the same conic section.

    Takes an optional `simplifier` callback that simplifies the incidence
    polynomial before it gets compared to zero. Returns `None` if undecidable.

    *Formula*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 10.2 (Conics and Cross-Ratios)
    """
    if len(points) < 6:
        return True
    if len(points) > 6:
        raise NotImplementedError
    points = [point_to_vec3(p) for p in points]
    d = []
    for p in points[0:2]:
        for q in points[2:4]:
            for r in points[4:6]:
                d.append(Matrix.hstack(p, q, r).det())  # noqa: PERF401
    incidence_poly = d[0] * d[3] * d[5] * d[6] - d[1] * d[2] * d[4] * d[7]
    return simplifier(incidence_poly).is_zero


def are_cocircular(
    points: Sequence[Matrix],
    *,
    simplifier: Callable[[Expr], Expr] = expand,
) -> bool | None:
    """Tells whether n points lie on the same circle.

    Note that four collinear points, as well as three collinear points and an
    arbitrary ideal point are also considered cocircular.

    Takes an optional `simplifier` callback that simplifies the cocircularity
    determinant before it gets compared to zero. Returns `None` if undecidable.

    *Formula*: [Measuring cocircularity in a point set](
    https://egc23.web.uah.es/wp-content/uploads/2023/06/EGC23_paper_20.pdf)
    *Own research*:
    [research/cocircularity.py](../src/research/cocircularity.py)
    """
    if len(points) < 4:
        return True

    cc_matrix_elements = []
    for p in points:
        x, y, z = point_to_vec3(p)
        cc_matrix_elements.append([x * z, y * z, x**2 + y**2, z**2])

    cocircularity_matrix = Matrix(cc_matrix_elements)
    if len(points) > 4:
        cocircularity_matrix = cocircularity_matrix.T * cocircularity_matrix

    return simplifier(cocircularity_matrix.det()).is_zero
