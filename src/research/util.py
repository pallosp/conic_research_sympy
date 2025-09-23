from sympy import Expr, pretty


def print_indented(expr: Expr, *, indent: int = 2) -> None:
    print(indent(pretty(expr), " " * indent))


def println_indented(expr: Expr, *, indent: int = 2) -> None:
    print(indent(pretty(expr), " " * indent) + "\n")
