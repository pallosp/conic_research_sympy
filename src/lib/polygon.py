from sympy import Expr, S

from lib.point import point_to_xy
from lib.type_aliases import FinitePoint


def signed_area(*points: FinitePoint) -> Expr:
    """Computes the signed area of the polygon with the given vertices.

    The area is positive if the vertices are in counterclockwise order, and
    negative if they are in clockwise order. For a self-intersecting polygon,
    it is the sum of the signed areas of its loops. It is zero for fewer than
    three points.

    *Formula*: <https://en.wikipedia.org/wiki/Shoelace_formula>
    """
    xy = [point_to_xy(p) for p in points]
    double_area = S.Zero
    for (x1, y1), (x2, y2) in zip(xy, xy[1:] + xy[:1], strict=True):
        double_area += x1 * y2 - x2 * y1
    return double_area / 2
