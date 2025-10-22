from collections.abc import Sequence

from sympy import Abs, Expr, I, Matrix, sqrt

from lib.conic_classification import ConicNormFactor
from lib.matrix import conic_matrix, max_eigenvalue, min_eigenvalue
from lib.point import point_to_xy


def conic_from_foci_and_radius(
    focus1: Matrix | Sequence[Expr],
    focus2: Matrix | Sequence[Expr],
    radius: Expr,
) -> Matrix:
    """Computes the ellipse or hyperbola with the given focus points and
    radius, i.e. center-vertex distance.

    If `radius` is negative, takes its absolute value.

    *Formula*:
    [research/conic_from_foci_and_radius.py](../src/research/conic_from_foci_and_radius.py)
    """
    fx1, fy1 = point_to_xy(focus1)
    fx2, fy2 = point_to_xy(focus2)

    # center
    cx, cy = (fx1 + fx2) / 2, (fy1 + fy2) / 2

    # center -> focus vector
    dx, dy = (fx2 - fx1) / 2, (fy2 - fy1) / 2

    a = dx**2 - radius**2
    b = dx * dy
    c = dy**2 - radius**2
    d = -cx * a - cy * b
    e = -cx * b - cy * c
    f = -cx * d - cy * e + (a * c - b * b)

    return conic_matrix(a, b, c, d, e, f)


def conic_from_center_and_points(
    center: Matrix | Sequence[Expr],
    p1: Matrix | Sequence[Expr],
    p2: Matrix | Sequence[Expr],
    p3: Matrix | Sequence[Expr],
) -> Matrix:
    """Computes the conic section with the given center and perimeter points.

    May return
     - an ellipse;
     - a hyperbola;
     - a parallel line pair;
     - zero matrix if the solution is ambiguous, which happens when some of the
       (`center`, `pᵢ`, `pⱼ`) triples are collinear.

    *Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)
    """
    x, y = point_to_xy(center)
    x1, y1 = point_to_xy(p1)
    x2, y2 = point_to_xy(p2)
    x3, y3 = point_to_xy(p3)
    x1, x2, x3 = x1 - x, x2 - x, x3 - x
    y1, y2, y3 = y1 - y, y2 - y, y3 - y
    m = Matrix(
        [
            [x1 * x1, x1 * y1, y1 * y1],
            [x2 * x2, x2 * y2, y2 * y2],
            [x3 * x3, x3 * y3, y3 * y3],
        ],
    )
    m1, m2, m3 = m.copy(), m.copy(), m.copy()
    m1[:, 0] = Matrix([1, 1, 1])
    m2[:, 1] = Matrix([1, 1, 1])
    m3[:, 2] = Matrix([1, 1, 1])
    a = m1.det()
    b = m2.det() / 2
    c = m3.det()
    d = -a * x - b * y
    e = -b * x - c * y
    f = -d * x - e * y - m.det()
    return conic_matrix(a, b, c, d, e, f)


def conic_center(conic: Matrix) -> Matrix:
    """Computes the center point of a conic.

    Returns the point's coordinates as a 2D column vector.

    *Formula*: [research/conic_center.py](../src/research/conic_center.py)
    """
    x, y, z = conic.col(0).cross(conic.col(1))
    return Matrix([x / z, y / z])


def semi_axis_lengths(conic: Matrix) -> tuple[Expr, Expr]:
    """Computes the semi-axis lengths of a conic in no specific order.

    To get the semi-focal or semi-transverse axis length (semi-major /
    semi-minor in case of ellipses), call
    [principal_radius](#central_conic.principal_radius) or
    [secondary_radius](#central_conic.secondary_radius), respectively.

    *Formula*: [research/conic_radii.py](../src/research/conic_radii.py)
    """
    submatrix = conic[:2, :2]
    return (
        sqrt(-conic.det() / (min_eigenvalue(submatrix) * submatrix.det())),
        sqrt(-conic.det() / (max_eigenvalue(submatrix) * submatrix.det())),
    )


def _selected_radius(conic: Matrix, eigenvalue_selector: Expr) -> Matrix:
    """Returns the conic radius corresponding to the selected eigenvalue.

    `eigenvalue_selector` is an expression evaluating to ±1. When positive, the
    larger, when negative the smaller conic eigenvalue will be used in the
    radius formula.
    """
    a, _, _, b, c, _, _, _, _ = conic
    eigenvalue = (a + c + eigenvalue_selector * sqrt((a - c) ** 2 + 4 * b**2)) / 2
    return sqrt(-conic.det() / (eigenvalue * (a * c - b * b)))


def primary_radius(conic: Matrix) -> Expr:
    """Computes the center-vertex distance of a conic.

    This corresponds to the semi-major axis length of real ellipses. In case of
    [imaginary ellipses](#conic_classification.is_imaginary_ellipse) however the
    focal axis is the shorter one in terms of absolute value.

    The returned value is:
     - a positive number for ellipses and hyperbolas;
     - infinity for parabolas;
     - an imaginary number for imaginary ellipses;
     - `nan` for ideal point conics;
     - 0 for the other degenerate conics.
    """
    return _selected_radius(conic, ConicNormFactor(conic))


def secondary_radius(conic: Matrix) -> Expr:
    """Computes the semi-conjugate axis length of a conic.

    This corresponds to the semi-minor axis length of real ellipses. In case of
    of imaginary ellipses however the conjugate axis is the longer one in terms
    absolute value. Hyperbolas intersect their conjugate axis at complex points,
    therefore the secondary radius will be a complex number.

    The returned value is:
     - a positive number for ellipses;
     - infinity for parabolas;
     - an imaginary number for hyperbolas and imaginary ellipses;
     - `nan` for ideal point conics;
     - 0 for the other degenerate conics.
    """
    return _selected_radius(conic, -ConicNormFactor(conic))


def linear_eccentricity(conic: Matrix) -> Expr:
    """Computes the linear eccentricity of a conic section.

    The linear eccentricity is the distance between the center and a focus
    point.

    Special cases:
     - zero for circles and imaginary circles;
     - a real number for imaginary ellipses, since they still have real center
       and foci;
     - infinity for parabolas;
     - `nan` for parallel and coincident line pairs;
     - `nan` for conics containing the ideal line;
     - zero for all other degenerate conics.

    *Formula*: `√|r₁²-r₂²|` where `r₁` and `r₂` denote the primary and secondary
    radii of the conic (i.e., the semi-axis lengths in the case of an ellipse).
    """
    a, _, _, b, c, _, _, _, _ = conic
    eigenvalue_diff = sqrt((a - c) ** 2 + 4 * b**2)
    return sqrt(Abs(conic.det()) * eigenvalue_diff) / Abs(a * c - b * b)


def center_to_focus_vector(conic: Matrix) -> Matrix:
    """Returns the 2D vector from a conic's center to one of its foci.

    The opposite vector points to the other focus.

    The function is only meaningful for
    [central conics](#conic_classification.is_central_conic). The result vector
    will contain infinite or `nan` elements when the conic lacks a finite
    center.
    """
    # Calculate the focal axis direction vector.
    a, _, _, b, c, _, _, _, _ = conic
    norm_sign = ConicNormFactor(conic)
    x, y = sqrt(norm_sign * (a - c + 2 * I * b)).simplify().as_real_imag()

    # Center-to-focus vector = [x, y] / √(x² + y²) * linear eccentricity
    # The √(x² + y²) = ∜((a-c)² + 4b²) factor vanishes.
    multiplier = sqrt(Abs(conic.det())) / (a * c - b * b)
    return Matrix([x, y]).applyfunc(lambda coord: coord * multiplier)


def center_to_vertex_vector(conic: Matrix) -> Matrix:
    """Returns the 2D vector from a conic's center to one of its vertices.

    The opposite vector points to the other vertex.

    The function is only meaningful for
    [central conics](#conic_classification.is_central_conic). The result vector
    will contain infinite or `nan` elements when the conic lacks a finite
    center.
    """
    # Calculate the focal axis direction vector and its length.
    a, _, _, b, c, _, _, _, _ = conic
    norm_sign = ConicNormFactor(conic)
    x, y = sqrt(norm_sign * (a - c + 2 * I * b)).simplify().as_real_imag()

    # Center-to-vertex vector = [x, y] / √(x² + y²) * primary radius
    # √(x² + y²) = ∜((a-c)² + 4b²)
    xy_length = sqrt(sqrt((a - c) ** 2 + 4 * b**2))
    multiplier = primary_radius(conic) / xy_length
    return Matrix([x, y]).applyfunc(lambda coord: coord * multiplier)


def shrink_conic_to_zero(conic: Matrix) -> Matrix:
    """Scales a conic section from its center with a factor of zero.

    Turns hyperbolas to line pair conics consisting of their asymptotes, and
    ellipses to point conics.

    This transformation is only meaningful for
    [central conics](#conic_classification.is_central_conic): for other
    conic types the result matrix will have infinite or `nan` elements.

    *Formula*:
    [research/scale_conic_from_center.py](../src/research/scale_conic_from_center.py
    """
    a, _, _, b, c, _, _, _, _ = conic
    return conic - Matrix.diag([0, 0, conic.det() / (a * c - b * b)])
