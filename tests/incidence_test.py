from sympy import symbols

from lib.incidence import AreCollinear


class TestAreCollinear:

    def test_less_than_three_points(self):
        assert AreCollinear((1, 2)) is True
        assert AreCollinear((1, 2), (3, 4)) is True
        assert AreCollinear((1, 2), (3, 4, 0)) is True

    def test_three_points_numeric(self):
        assert AreCollinear((1, 2), (3, 4), (5, 6)) is True
        assert AreCollinear((1, 2), (3, 4), (5, 7)) is False
        assert AreCollinear((1, 2), (3, 4), (1, 1, 0)) is True
        assert AreCollinear((1, 2), (3, 4), (1, 2, 0)) is False

    def test_three_points_symbolic(self):
        x = symbols("x")
        assert AreCollinear((1, 2), (3, 4), (x, x + 1)) is True
        assert AreCollinear((1, 2), (3, 4), (x, x)) is False
        assert AreCollinear((1, 2), (3, 4), symbols("x y")) is None

    def test_four_points_numeric(self):
        assert AreCollinear((1, 2), (1, 2), (3, 4), (5, 6)) is True
        assert AreCollinear((1, 2), (1, 2), (3, 4), (5, 7)) is False
        assert AreCollinear((1, 2), (3, 4), (5, 6), (7, 8)) is True
        assert AreCollinear((1, 2), (3, 4), (5, 6), (7, 9)) is False
