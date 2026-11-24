import pytest
from sympy import Expr, I, Piecewise, S, Symbol, im, re

from lib.math import sqrt_complex


class TestSqrtComplex:
    # (z, sqrt(z)) pairs covering zero, both real half-axes and both half-planes.
    SQUARE_ROOTS = (
        (S(0), S(0)),
        (S(4), S(2)),
        (S(-4), 2 * I),
        (3 + 4 * I, 2 + I),
        (3 - 4 * I, 2 - I),
    )

    @pytest.mark.parametrize(("z", "root"), SQUARE_ROOTS)
    def test_numbers(self, z: Expr, root: Expr):
        x, y = sqrt_complex(z)
        assert x == re(root)
        assert y == im(root)

    @pytest.mark.parametrize(("z_value", "root"), SQUARE_ROOTS)
    def test_substituted_later(self, z_value: Expr, root: Expr):
        z = Symbol("z")
        x, y = sqrt_complex(z)
        assert y.has(Piecewise)
        assert x.subs(z, z_value) == re(root)
        assert y.subs(z, z_value) == im(root)
