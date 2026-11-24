from sympy import Expr, Piecewise, sqrt


def sqrt_complex(z: Expr) -> tuple[Expr, Expr]:
    """Computes the principal square root of a complex number.

    Returns `(re(√z), im(√z))` expressed with `sqrt` and `Piecewise`.
    If the sign of `im(z)` is unknown, `Piecewise` stays unevaluated.

    Unlike `sqrt(z).as_real_imag()`, which goes through `atan2(im(z), re(z))`,
    it doesn't become `nan` when `z` is substituted with 0 later.
    See https://github.com/sympy/sympy/issues/29000 for more details.
    """
    x, y = z.as_real_imag()
    abs_z = sqrt(x * x + y * y)
    im_sign = Piecewise((1, y >= 0), (-1, True))
    return (sqrt((abs_z + x) / 2), im_sign * sqrt((abs_z - x) / 2))
