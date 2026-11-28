from textwrap import indent as indent_multiline

from sympy import Basic
from sympy.printing.pretty.pretty import PrettyPrinter

HORIZONTAL_LINE = "=" * 88


class _ResearchPrettyPrinter(PrettyPrinter):
    """Pretty printer with a longer implication arrow."""

    def _print_Implies(self, e: Basic, altchar: str | None = None) -> object:  # noqa: N802
        return super()._print_Implies(e, altchar or " ─→ ")


def _pretty(expr: object) -> str:
    return _ResearchPrettyPrinter().doprint(expr)


def print_indented(expr: object, *, indent: int = 2) -> None:
    print(indent_multiline(_pretty(expr), " " * indent))


def println_indented(expr: object, *, indent: int = 2) -> None:
    print(indent_multiline(_pretty(expr), " " * indent) + "\n")
