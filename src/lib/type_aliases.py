from typing import TypeAlias

from sympy import Expr

#: A scalar value: a SymPy expression or a Python number. Functions convert
#: numbers to SymPy expressions as needed.
Scalar: TypeAlias = Expr | int | float
