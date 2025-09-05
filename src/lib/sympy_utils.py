from sympy import Eq, Expr, Pow


def AddEq(*eqs: Eq) -> Eq:
    """Adds multiple sympy equations."""
    left = eqs[0].lhs
    right = eqs[0].rhs
    for term in eqs[1:]:
        left += term.lhs
        right += term.rhs
    return Eq(left, right)


def SubEq(eq0: Eq, eq1: Eq) -> Eq:
    """Subtracts one sympy equation from another."""
    return Eq(eq0.lhs - eq1.lhs, eq0.rhs - eq1.rhs)


def MulEq(eq: Eq, factor: Expr) -> Eq:
    """Multiplies a sympy equation by a factor."""
    return Eq(eq.lhs * factor, eq.rhs * factor)


def DivEq(eq: Eq, denom: Expr) -> Eq:
    """Divides a sympy equation by a denominator."""
    return Eq(eq.lhs / denom, eq.rhs / denom)


def SwapEq(eq: Eq) -> Eq:
    """Swaps the lhs and rhs of a sympy equation."""
    return Eq(eq.rhs, eq.lhs)


def FactorRadicals(expr: Expr) -> Expr:
    """Factors all `Pow` and `sqrt` subexpressions inside `expr`."""
    return expr.replace(Pow, lambda base, exp: Pow(base.factor(), exp))
