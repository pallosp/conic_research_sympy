import pytest
from sympy import Determinant, Matrix, S, sqrt, symbols, zeros

from lib.polygon import signed_area
from research.poly_as_det import poly_as_det3, poly_as_det3_matrix

x0, x1, x2, x3, x5 = symbols("x0 x1 x2 x3 x5")
y0, y1, y2, y3, y5 = symbols("y0 y1 y2 y3 y5")
u0, u1, u2, u3 = symbols("u0 u1 u2 u3")
v1, v2, v3 = symbols("v1 v2 v3")
w1, w2, w3 = symbols("w1 w2 w3")
c1, c4, c5, c6 = symbols("c1 c4 c5 c6")


class TestPolyAsDet3:
    def test_consecutive_indices(self):
        poly = signed_area((x1, y1), (x2, y2), (x3, y3))
        det = poly_as_det3(poly)
        assert det == Determinant(Matrix([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]])) / 2

    def test_non_consecutive_indices(self):
        poly = 2 * signed_area((x0, y0), (x2, y2), (x5, y5))
        assert poly_as_det3(poly) == Determinant(
            Matrix([[x0, y0, 1], [x2, y2, 1], [x5, y5, 1]])
        )

    def test_factor_with_other_index(self):
        det = Matrix([[u1, v1, w1], [u2, v2, w2], [u3, v3, w3]]).det()
        k, matrix = poly_as_det3_matrix(u0 * det)
        assert (k * matrix.det()).expand() == (u0 * det).expand()
        assert k in (u0, -u0)
        assert u0 not in matrix.free_symbols

    def test_factor_without_other_indices(self):
        poly = 2 * c1 * signed_area((x1, y1), (x2, y2), (x3, y3))
        k, matrix = poly_as_det3_matrix(poly)
        assert (k * matrix.det()).expand() == poly.expand()
        assert k in (c1, -c1)
        assert c1 not in matrix.free_symbols

    def test_several_index_shifts(self):
        # Invariant under 1 → 2 → 3 → 1 and 4 → 5 → 6 → 4
        poly = 2 * (c4 + c5 + c6) * signed_area((x1, y1), (x2, y2), (x3, y3))
        k, matrix = poly_as_det3_matrix(poly)
        assert (k * matrix.det()).expand() == poly.expand()
        assert k == 1  # the non-monomial factor goes into the matrix

    def test_zero(self):
        assert poly_as_det3_matrix(S.Zero) == (0, zeros(3, 3))


class TestPolyAsDet3Errors:
    def test_irrational_coefficient(self):
        with pytest.raises(ValueError, match="rational coefficients"):
            poly_as_det3(2 * sqrt(2) * signed_area((x1, y1), (x2, y2), (x3, y3)))

    def test_monomial(self):
        with pytest.raises(ValueError, match="cannot express"):
            poly_as_det3(6 * x1 * x2 * x3)

    def test_not_cyclic(self):
        with pytest.raises(ValueError, match="not invariant under any"):
            poly_as_det3(x1 * y2 + x2 * y3)

    def test_several_index_shifts_without_determinant(self):
        with pytest.raises(ValueError, match="cannot express"):
            poly_as_det3(x0 + x1 + x2 + x3)

    def test_constants_hint(self):
        with pytest.raises(ValueError, match=r"treated as constants: c1\)"):
            poly_as_det3(u1 * (x1 + x2 + x3 + c1))
