from textwrap import indent as indent_multiline

from sympy import Expr, pretty


def print_indented(expr: Expr, *, indent: int = 2) -> None:
    print(indent_multiline(pretty(expr), " " * indent))


def println_indented(expr: Expr, *, indent: int = 2) -> None:
    print(indent_multiline(pretty(expr), " " * indent) + "\n")
