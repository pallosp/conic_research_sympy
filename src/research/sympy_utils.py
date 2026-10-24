from textwrap import indent as indent_multiline
from typing import TYPE_CHECKING, TypeAlias, cast

from sympy import Basic, Eq, Equality, Expr, MatrixBase, Pow, pretty

from lib.type_aliases import Scalar

if TYPE_CHECKING:

    class ExprEq(Equality):
        """An equation whose sides are expressions.

        SymPy types the sides of `Equality` as `Basic`, which doesn't support
        arithmetic. Exists only for type checking; at runtime it's `Equality`.
        """

        @property
        def lhs(self) -> Expr:
            """The left-hand side of the equation."""
            ...

        @property
        def rhs(self) -> Expr:
            """The right-hand side of the equation."""
            ...

else:
    ExprEq = Equality


def eq(lhs: Scalar, rhs: Scalar) -> ExprEq:
    """Creates a sympy equation without evaluating it to True or False."""
    return cast("ExprEq", Eq(lhs, rhs, evaluate=False))


def add_eq(*eqs: ExprEq) -> ExprEq:
    """Adds multiple sympy equations."""
    left = eqs[0].lhs
    right = eqs[0].rhs
    for term in eqs[1:]:
        left += term.lhs
        right += term.rhs
    return eq(left, right)


def sub_eq(eq0: ExprEq, eq1: Expr | ExprEq) -> ExprEq:
    """Subtracts a sympy equation or expression from another equation."""
    lhs = eq1.lhs if isinstance(eq1, ExprEq) else eq1
    rhs = eq1.rhs if isinstance(eq1, ExprEq) else eq1
    return eq(eq0.lhs - lhs, eq0.rhs - rhs)


def mul_eq(equation: ExprEq, factor: Scalar | ExprEq) -> ExprEq:
    """Multiplies a sympy equation by a factor or another equation."""
    lhs = factor.lhs if isinstance(factor, ExprEq) else factor
    rhs = factor.rhs if isinstance(factor, ExprEq) else factor
    return eq(equation.lhs * lhs, equation.rhs * rhs)


def div_eq(equation: ExprEq, denom: Scalar) -> ExprEq:
    """Divides a sympy equation by a denominator."""
    return eq(equation.lhs / denom, equation.rhs / denom)


def swap_eq(equation: ExprEq) -> ExprEq:
    """Swaps the lhs and rhs of a sympy equation."""
    return eq(equation.rhs, equation.lhs)


#: An item of an equation chain built by `eq_chain`.
EqChainItem: TypeAlias = Basic | MatrixBase | tuple[object, ...] | int | float


def eq_chain(*expressions: EqChainItem) -> Eq | EqChainItem:
    """Creates a chain of equations, i.e. `expr_1 = expr_2 = ... = expr_n`."""
    if not expressions:
        raise ValueError("At least one expression is required")
    if len(expressions) == 1:
        return expressions[0]
    chain = Eq(eq_chain(*expressions[0:-1]), expressions[-1], evaluate=False)
    return cast("Eq", chain)


def factor_radicals(expr: Expr) -> Expr:
    """Factors all `Pow` and `sqrt` subexpressions inside `expr`.

    Note that there is a simpler alternative that factors all subexpressions:
    `expr.factor(deep=True)`
    """
    return expr.subs({p: Pow(p.base.factor(), p.exp) for p in expr.atoms(Pow)})


def print_indented(expr: object, *, indent: int = 2) -> None:
    print(indent_multiline(pretty(expr), " " * indent))


def println_indented(expr: object, *, indent: int = 2) -> None:
    print(indent_multiline(pretty(expr), " " * indent) + "\n")
