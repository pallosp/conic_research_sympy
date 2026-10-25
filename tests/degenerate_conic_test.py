from typing import TYPE_CHECKING, cast

import pytest
from sympy import I, Matrix, MatrixBase, symbols

from lib.circle import circle
from lib.conic import conic_from_poly, ideal_points
from lib.conic_classes import is_degenerate, is_finite_conic
from lib.degenerate_conic import (
    ZERO_CONIC,
    double_line_conic,
    extract_point,
    line_pair_conic,
    point_conic,
    split_to_lines,
)
from lib.line import (
    IDEAL_LINE,
    X_AXIS,
    Y_AXIS,
    ZERO_LINE,
    horizontal_line,
    vertical_line,
)
from lib.matrix import conic_matrix, is_nonzero_multiple, is_real_matrix
from lib.point import ZERO_POINT, point_to_vec3
from tests.utils import are_projective_sets_equal, is_explicit_matrix

if TYPE_CHECKING:
    from lib.type_aliases import Line


class TestLinePair:
    def test_invalid_argument(self):
        with pytest.raises(ValueError, match="must be 3-dimensional column vectors"):
            line_pair_conic(Matrix([1, 2, 3]), Matrix([1, 2, 3]).T)
        with pytest.raises(ValueError, match="must be 3-dimensional column vectors"):
            line_pair_conic(Matrix([1, 2, 3]).T, Matrix([1, 2, 3]))

    def test_numeric_line_pair(self):
        double_ideal = conic_matrix(0, 0, 0, 0, 0, 1)
        assert is_nonzero_multiple(
            line_pair_conic(IDEAL_LINE, IDEAL_LINE),
            double_ideal,
        )
        plus = conic_matrix(0, 1, 0, 0, 0, 0)
        assert is_nonzero_multiple(line_pair_conic(X_AXIS, Y_AXIS), plus)


class TestDoubleLine:
    def test_invalid_argument(self):
        with pytest.raises(ValueError, match="must be a 3-dimensional column vector"):
            double_line_conic(Matrix([1, 2, 3]).T)

    def test_numeric_line_pair(self):
        assert double_line_conic(X_AXIS) == conic_from_poly(symbols("y") ** 2)


class TestPointConic:
    def test_real_point(self):
        point = (1, 2)
        conic = point_conic(point)
        assert is_degenerate(conic) is True
        assert is_finite_conic(conic) is True

    def test_ideal_point(self):
        point = (1, 2, 0)
        conic = point_conic(point)
        assert is_degenerate(conic) is True
        assert are_projective_sets_equal(ideal_points(conic), [point, point])
        line = split_to_lines(conic)[0]
        assert is_explicit_matrix(line)
        assert is_real_matrix(line) is False


class TestSplitToLines:
    @pytest.mark.parametrize(
        "lines",
        [
            [horizontal_line(1), horizontal_line(2)],
            [vertical_line(1), vertical_line(2)],
            [Matrix([1, 2, 3])] * 2,
            [IDEAL_LINE] * 2,
            [Matrix([1, 2, 3]), Matrix([4, 5, 6])],
        ],
        ids=["horizontal", "vertical", "double", "double_ideal", "general"],
    )
    def test_line_pair(self, lines: list[MatrixBase]):
        conic = line_pair_conic(*lines)
        assert are_projective_sets_equal(lines, split_to_lines(conic))

    def test_point_conic(self):
        conic = circle((0, 0), r=0)
        assert are_projective_sets_equal(
            [Matrix([1, I, 0]), Matrix([1, -I, 0])],
            split_to_lines(conic),
        )

    def test_zero_conic(self):
        assert split_to_lines(ZERO_CONIC) == (ZERO_LINE, ZERO_LINE)

    def test_symbolic_conic(self):
        symbolic_conic = conic_matrix(*symbols("a b c d e f"))
        lines = split_to_lines(symbolic_conic)
        assert not is_explicit_matrix(lines[0])
        assert not is_explicit_matrix(lines[1])

        # Substituting a concrete line pair conic evaluates them to its lines.
        line_pair = [Matrix([1, 2, 3]), Matrix([4, 5, 6])]
        conic = line_pair_conic(*line_pair)
        substitutions = list(zip(symbolic_conic.flat(), conic.flat(), strict=True))
        concrete_lines = [line.subs(substitutions) for line in lines]
        assert is_explicit_matrix(concrete_lines[0])
        assert is_explicit_matrix(concrete_lines[1])
        # The lines are explicit, as asserted above. SymPy types MatrixExpr.subs()
        # as returning Expr (fix upstream).
        assert are_projective_sets_equal(cast("list[Line]", concrete_lines), line_pair)

    def test_partially_symbolic_conic(self):
        symbolic_conic = conic_matrix(1, *symbols("b c d e f"))
        lines = split_to_lines(symbolic_conic)
        # The intersection point is undecidable, but the pivot isn't: the lines
        # are explicit matrices with unevaluated elements.
        assert is_explicit_matrix(lines[0])
        assert is_explicit_matrix(lines[1])

        # Substituting a concrete line pair conic evaluates them to its lines.
        # Its x² coefficient is 1·1 = 1, matching the symbolic conic.
        line_pair = [Matrix([1, 2, 3]), Matrix([1, 5, 6])]
        conic = line_pair_conic(*line_pair)
        substitutions = list(zip(symbolic_conic.flat(), conic.flat(), strict=True))
        concrete_lines = [line.subs(substitutions) for line in lines]
        assert is_explicit_matrix(concrete_lines[0])
        assert is_explicit_matrix(concrete_lines[1])
        # The lines are explicit, as asserted above. SymPy types MatrixExpr.subs()
        # as returning Expr (fix upstream).
        assert are_projective_sets_equal(cast("list[Line]", concrete_lines), line_pair)


class TestExtractPoint:
    def test_symbolic_real_point_conic(self):
        point = symbols("x,y", real=True)
        zero_circle = circle(point, r=0)
        assert is_nonzero_multiple(extract_point(zero_circle), point_to_vec3(point))

        finite_point_conic = point_conic(point)
        assert is_nonzero_multiple(
            extract_point(finite_point_conic),
            point_to_vec3(point),
        )

    def test_symbolic_ideal_point_conic(self):
        point = [*symbols("x,y", positive=True), 0]
        conic = point_conic(point)
        assert is_nonzero_multiple(extract_point(conic), point_to_vec3(point))

    def test_zero_matrix(self):
        assert extract_point(ZERO_CONIC) == ZERO_POINT

    def test_numeric_line_pairs(self):
        line1 = Matrix([1, 2, 3])
        line2 = Matrix([4, 5, 6])
        conic = line_pair_conic(line1, line2)
        intersection = line1.cross(line2)
        assert is_nonzero_multiple(extract_point(conic), intersection)

    def test_symbolic_double_line(self):
        line = Matrix(symbols("a b c"))
        conic = line_pair_conic(line, line)
        assert extract_point(conic) == ZERO_POINT
