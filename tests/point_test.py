import pytest
from sympy import I, Matrix, MatrixBase, MatrixSymbol, Rational, nan, symbols, zoo

from lib.circle import circle
from lib.conic_classes import is_circular
from lib.degenerate_conic import double_line_conic, point_conic
from lib.incidence import conic_contains_point
from lib.line import IDEAL_LINE, X_AXIS
from lib.matrix import is_nonzero_multiple, pivot_column
from lib.point import (
    CIRCULAR_POINT_I,
    CIRCULAR_POINT_J,
    ZERO_POINT,
    centroid,
    direction_to_xy,
    ideal_point_on_line,
    perpendicular_foot,
    point_to_vec3,
    point_to_xy,
)
from tests.utils import is_explicit_matrix


class TestPointToVec3:
    def test_2d(self):
        assert point_to_vec3((1, 2)) == Matrix([1, 2, 1])
        assert point_to_vec3(Matrix([1, 2])) == Matrix([1, 2, 1])

    def test_3d(self):
        assert point_to_vec3((1, 2, 3)) == Matrix([1, 2, 3])
        assert point_to_vec3([1, 2, 3]) == Matrix([1, 2, 3])

    def test_column_vector_is_not_copied(self):
        point = Matrix([1, 2, 3])
        assert point_to_vec3(point) is point

    def test_row_vector(self):
        assert point_to_vec3(Matrix([[1, 2, 3]])) == Matrix([1, 2, 3])

    def test_invalid_length(self):
        with pytest.raises(ValueError, match="2D or 3D"):
            point_to_vec3([1, 2, 3, 4])

    def test_invalid_shape(self):
        with pytest.raises(ValueError, match="2D or 3D"):
            point_to_vec3(Matrix([[1, 2], [3, 4]]))

    def test_matrix_expr(self):
        point = MatrixSymbol("p", 3, 1)
        assert point_to_vec3(point) is point
        row = MatrixSymbol("r", 1, 3)
        assert point_to_vec3(row) == row.T
        xy = MatrixSymbol("q", 2, 1)
        assert point_to_vec3(xy) == Matrix([xy[0, 0], xy[1, 0], 1])

    def test_pivot_column(self):
        column = pivot_column(Matrix(3, 3, symbols("m:9")))
        assert not is_explicit_matrix(column)
        assert point_to_vec3(column) is column


class TestPointToXy:
    def test_2d(self):
        assert point_to_xy((1, 2)) == (1, 2)
        assert point_to_xy([1, 2]) == (1, 2)
        assert point_to_xy(Matrix([1, 2])) == (1, 2)

    def test_3d(self):
        assert point_to_xy((1, 2, 3)) == (Rational(1, 3), Rational(2, 3))
        assert point_to_xy([-3, -2, -1]) == (3, 2)
        assert point_to_xy(Matrix([-3, -2, -1])) == (3, 2)

    def test_ideal_point(self):
        assert point_to_xy((1, 2, 0)) == (zoo, zoo)

    def test_zero_vector(self):
        assert point_to_xy((0, 0, 0)) == (nan, nan)


class TestDirectionToXy:
    def test_2d(self):
        assert direction_to_xy((1, 2)) == (1, 2)
        assert direction_to_xy(Matrix([1, 2])) == (1, 2)

    def test_ideal_point(self):
        assert direction_to_xy((1, 2, 0)) == (1, 2)
        assert direction_to_xy(Matrix([1, 2, 0])) == (1, 2)

    def test_finite_point(self):
        with pytest.raises(ValueError, match="must be a 2D vector or an ideal point"):
            direction_to_xy((1, 2, 1))

    def test_wrong_length(self):
        with pytest.raises(ValueError, match="must be a 2D vector or an ideal point"):
            direction_to_xy([1, 2, 0, 0])


class TestIdealPointOnLine:
    def test_symbolic(self):
        line = Matrix(symbols("a b c"))
        point = ideal_point_on_line(line)
        assert point[2] == 0  # point should be an ideal point
        assert line.dot(point) == 0  # point should be on the line

    def test_ideal_line(self):
        assert ideal_point_on_line(IDEAL_LINE) == ZERO_POINT


class TestCentroid:
    def test_no_points(self):
        with pytest.raises(ValueError, match="At least one point"):
            centroid()

    def test_one_point(self):
        assert centroid((1, 2)) == Matrix([1, 2])
        assert centroid((6, 4, 2)) == Matrix([3, 2])

    def test_two_points(self):
        assert centroid((1, 2), (3, 4)) == Matrix([2, 3])

    def test_full_precision(self):
        expected = Matrix([Rational(2, 3), Rational(1, 3)])
        assert centroid((0, 0), (1, 0), (1, 1)) == expected


class TestPerpendicularFoot:
    def test_finite_point_and_line(self):
        point = (1, -1)
        line = Matrix([1, 1, -4])  # x+y=4
        assert is_nonzero_multiple(perpendicular_foot(point, line), (3, 1))

    def test_ideal_point_or_line(self):
        nan_point = Matrix([nan, nan])
        assert perpendicular_foot((1, 2), IDEAL_LINE) == nan_point
        assert perpendicular_foot((1, 2, 0), Matrix([1, 2, 3])) == nan_point
        assert perpendicular_foot((1, 2, 0), IDEAL_LINE) == nan_point


class TestCircularPoints:
    @pytest.mark.parametrize("point", [CIRCULAR_POINT_I, CIRCULAR_POINT_J])
    def test_on_circles(self, point: MatrixBase) -> None:
        assert conic_contains_point(circle((3, -2), r=5), point)
        assert conic_contains_point(circle((3, -2), r=0), point)
        assert conic_contains_point(circle((3, -2), r=2 * I), point)

    @pytest.mark.parametrize("point", [CIRCULAR_POINT_I, CIRCULAR_POINT_J])
    def test_on_double_ideal_line(self, point: MatrixBase) -> None:
        assert conic_contains_point(double_line_conic(IDEAL_LINE), point)

    @pytest.mark.parametrize("point", [CIRCULAR_POINT_I, CIRCULAR_POINT_J])
    def test_not_on_non_circular_point_conic(self, point: MatrixBase) -> None:
        conic = point_conic((3, -2))
        assert is_circular(conic) is False
        assert conic_contains_point(conic, point) is False

    @pytest.mark.parametrize("point", [CIRCULAR_POINT_I, CIRCULAR_POINT_J])
    def test_not_on_double_real_line(self, point: MatrixBase) -> None:
        assert conic_contains_point(double_line_conic(X_AXIS), point) is False
