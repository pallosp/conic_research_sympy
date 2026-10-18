from sympy import MatrixBase, MatrixExpr, sqrt

from lib.matrix import pivot_column_and_row, quadratic_form, skew_matrix
from lib.type_aliases import Line, PointVector


def line_x_line(line1: Line, line2: Line) -> PointVector:
    """Computes the intersection of two lines.

    Returns an ideal point if the lines are parallel, or `[0, 0, 0]ᵀ` if they
    coincide.
    """
    return line1.cross(line2)


def conic_x_line(
    conic: MatrixBase,
    line: Line,
) -> tuple[PointVector | MatrixExpr, PointVector | MatrixExpr]:
    """Intersects a conic with a line. Returns two points.

    Special cases:
     - The intersection points coincide if the line is tangent to the conic.
     - They are complex conjugates if the line doesn't intersect the conic at
       a real point.
     - May return unevaluated `MatrixExpr`s for symbolic conics or lines,
       when none of the intersection matrix elements can be proven non-zero.
     - Returns two `[0, 0, 0]ᵀ` vectors if the conic contains the entire line.

    *Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
    section 11.3
    """
    # Notation: l = line, C = conic, B = skew_mat.
    #
    # The book's algorithm:
    #  1. m = Bᵀ·C·B is a degenerate conic representing the two intersection
    #     points p and q as a symmetric matrix: m ~ p·qᵀ + q·pᵀ.
    #  2. Adding α·B with the right α splits it into the rank-1 matrix p·qᵀ.
    #     The book computes α from any non-zero coordinate of the line:
    #     α = √(-adj(m)ᵢᵢ) / lᵢ, where adj(m)ᵢᵢ is the determinant of the 2×2
    #     matrix left after deleting row i and column i of m.
    #  3. A non-zero element of p·qᵀ determines the points: its column is p,
    #     its row is q.
    #
    # Branch reduction: step 2 requires picking an i with lᵢ ≠ 0. We avoid it
    # by computing α without division:
    #  - adj(m) = adj(B)·adj(C)·adj(Bᵀ), since adj(X·Y) = adj(Y)·adj(X).
    #  - adj(B) = l·lᵀ for the skew matrix of l.
    #  - Hence adj(m) = (lᵀ·adj(C)·l)·l·lᵀ, so adj(m)ᵢᵢ / lᵢ² = lᵀ·adj(C)·l.
    #  - This gives α = √(-lᵀ·adj(C)·l) for any line.
    #
    # Only the branch in step 3 (pivot_column_and_row) remains, and it's
    # unavoidable.

    skew_mat = skew_matrix(line)
    m = skew_mat.T * conic * skew_mat
    alpha = sqrt(-quadratic_form(conic.adjugate(), line))
    intersections = m + alpha * skew_mat
    return pivot_column_and_row(intersections)
