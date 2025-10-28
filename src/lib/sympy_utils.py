from sympy import Eq, Expr


def add_eq(*eqs: Eq) -> Eq:
    """Adds multiple sympy equations."""
    left = eqs[0].lhs
    right = eqs[0].rhs
    for term in eqs[1:]:
        left += term.lhs
        right += term.rhs
    return Eq(left, right)


def sub_eq(eq0: Eq, eq1: Expr | Eq) -> Eq:
    """Subtracts a sympy equation or expression from another equation."""
    lhs = eq1.lhs if isinstance(eq1, Eq) else eq1
    rhs = eq1.rhs if isinstance(eq1, Eq) else eq1
    return Eq(eq0.lhs - lhs, eq0.rhs - rhs)


def mul_eq(eq: Eq, factor: Expr | Eq) -> Eq:
    """Multiplies a sympy equation by a factor or another equation."""
    lhs = factor.lhs if isinstance(factor, Eq) else factor
    rhs = factor.rhs if isinstance(factor, Eq) else factor
    return Eq(eq.lhs * lhs, eq.rhs * rhs)


def div_eq(eq: Eq, denom: Expr) -> Eq:
    """Divides a sympy equation by a denominator."""
    return Eq(eq.lhs / denom, eq.rhs / denom)


def swap_eq(eq: Eq) -> Eq:
    """Swaps the lhs and rhs of a sympy equation."""
    return Eq(eq.rhs, eq.lhs)


def eq_chain(*expressions: Expr) -> Eq | Expr:
    """Creates a chain of equations, i.e. `expr_1 = expr_2 = ... = expr_n`."""
    if not expressions:
        raise ValueError("At least one expression is required")
    if len(expressions) == 1:
        return expressions[0]
    return Eq(eq_chain(*expressions[0:-1]), expressions[-1], evaluate=False)
