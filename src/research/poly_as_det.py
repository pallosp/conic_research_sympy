"""Expresses a cyclic polynomial as a 3x3 determinant.

The rows of the matrix are cyclic shifts of each other: the next row is the
previous one with the variable indices shifted by 1 → 2 → 3 → 1. The entries
of a column have the same degree, non-increasing from left to right.

Search: start from a monomial matrix whose determinant is the closest to the
polynomial, then keep adding ±monomials to the entries of the first row as
long as that brings the determinant closer. The distance of two polynomials is
the sum of the absolute values of the coefficients of their difference. Ties
are broken randomly, and failed runs are retried.

Example:
    >>> from sympy import symbols
    >>> from research.poly_as_det import poly_as_det3
    >>> x1, x2, x3, y1, y2, y3 = symbols("x1 x2 x3 y1 y2 y3")
    >>> poly_as_det3(2 * (x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y1 - x1*y3))
    2*Determinant(Matrix([
    [x1, y1, 1],
    [x2, y2, 1],
    [x3, y3, 1]]))

"""

import itertools
import math
import random
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from functools import cached_property
from typing import TypeVar

from sympy import (
    ZZ,
    Determinant,
    Expr,
    Matrix,
    Poly,
    Rational,
    Symbol,
    sympify,
    zeros,
)
from sympy.combinatorics import Permutation
from sympy.polys.rings import PolyElement, PolyRing

Exponents = tuple[int, ...]  # exponents of a monomial
Row = tuple[PolyElement, PolyElement, PolyElement]  # first row of the matrix
T = TypeVar("T")


def _index_shift(symbols: Iterable[Symbol]) -> dict[Symbol, Symbol]:
    """Maps x1 → x2, x2 → x3, x3 → x1 for every prefix x with all 3 indices.

    Other symbols, such as k1 without k2 and k3, or x12, are mapped to
    themselves.
    """
    by_name = {s.name: s for s in symbols}
    shift = {s: s for s in by_name.values()}
    for name, symbol in by_name.items():
        match = re.fullmatch(r"(.*\D)([123])", name)
        if not match:
            continue
        prefix, index = match.groups()
        if all(f"{prefix}{i}" in by_name for i in (1, 2, 3)):
            shift[symbol] = by_name[f"{prefix}{int(index) % 3 + 1}"]
    return shift


def _divisors(exps: Exponents) -> Iterator[Exponents]:
    """Yields the exponents of all monomials dividing the given one."""
    return itertools.product(*(range(e + 1) for e in exps))


def _sub(a: Exponents, b: Exponents) -> Exponents:
    """Exponents of monomial a divided by monomial b."""
    return tuple(x - y for x, y in zip(a, b, strict=True))


def _distance(p: PolyElement) -> int:
    """Sum of the absolute values of the coefficients."""
    return sum(abs(c) for c in p.coeffs())


def _closest(scored: Iterable[tuple[int, T]], baseline: int) -> list[T]:
    """Items with the smallest distance, provided that it's below baseline."""
    best, items = baseline, []
    for distance, item in scored:
        if distance < best:
            best, items = distance, [item]
        elif distance == best and items:
            items.append(item)
    return items


class _Search:
    """Randomized greedy search for a cyclic matrix with a given determinant."""

    def __init__(
        self,
        target: PolyElement,
        index_shift: dict[Symbol, Symbol],
        max_degree: int | None,
    ) -> None:
        self.target = target
        self.ring = target.ring
        gens = list(self.ring.symbols)
        self.shift_perm = [gens.index(index_shift[g]) for g in gens]
        self.max_degree = max_degree
        # Candidate monomials for the matrix entries, grouped by degree
        self.monomials: dict[int, list[Exponents]] = defaultdict(list)
        for exps in sorted(
            {d for m in target.monoms() for d in _divisors(m) if self._fits(d)}
        ):
            self.monomials[sum(exps)].append(exps)

    def _fits(self, exps: Exponents) -> bool:
        """Tells whether the monomial is within the degree limit."""
        return self.max_degree is None or sum(exps) <= self.max_degree

    def _shift_exps(self, exps: Exponents, times: int = 1) -> Exponents:
        """Shifts the indices of a monomial's variables; times may be negative."""
        for _ in range(times % 3):
            shifted = [0] * len(exps)
            for i, e in enumerate(exps):
                shifted[self.shift_perm[i]] = e
            exps = tuple(shifted)
        return exps

    def shift(self, p: PolyElement, times: int = 1) -> PolyElement:
        """Shifts the indices of the variables in p: x1 → x2 → x3 → x1."""
        return self.ring.from_dict(
            {self._shift_exps(m, times): c for m, c in p.terms()}
        )

    def monomial(self, exps: Exponents, coeff: int = 1) -> PolyElement:
        return self.ring.from_dict({exps: coeff})

    def rows(self, row: Row) -> list[list[PolyElement]]:
        """The rows of the matrix: row, shift(row), shift²(row)."""
        return [[self.shift(e, i) for e in row] for i in range(3)]

    def cofactors(self, row: Row, col: int) -> list[PolyElement]:
        """Cofactors of the entries in column col of the matrix."""
        rows = self.rows(row)
        left, right = (j for j in range(3) if j != col)
        cofactors = []
        for i in range(3):
            upper, lower = (rows[k] for k in range(3) if k != i)
            minor = upper[left] * lower[right] - upper[right] * lower[left]
            cofactors.append(minor if (i + col) % 2 == 0 else -minor)
        return cofactors

    def det(self, row: Row) -> PolyElement:
        """Determinant of the matrix, expanded along the first column."""
        return sum(
            (self.shift(row[0], i) * c for i, c in enumerate(self.cofactors(row, 0))),
            self.ring.zero,
        )

    def _start_candidates(self) -> Iterator[Row]:
        """Yields the monomial rows (a, b, c) worth trying as a start.

        A useful start's diagonal term a⋅shift(b)⋅shift²(c) is in the target.
        The degrees are non-increasing: a column swap only flips the sign.
        """
        for term in self.target.monoms():
            for a in _divisors(term):
                if not self._fits(a):
                    continue
                for shifted_b in _divisors(_sub(term, a)):
                    shifted_c = _sub(_sub(term, a), shifted_b)
                    b = self._shift_exps(shifted_b, -1)
                    c = self._shift_exps(shifted_c, -2)
                    if sum(a) >= sum(b) >= sum(c):
                        for sign in (1, -1):
                            yield (
                                self.monomial(a, sign),
                                self.monomial(b),
                                self.monomial(c),
                            )

    @cached_property
    def starts(self) -> list[Row]:
        """Monomial rows whose determinant is the closest to the target."""
        scored = (
            (_distance(self.target - self.det(row)), row)
            for row in set(self._start_candidates())
        )
        return sorted(_closest(scored, _distance(self.target)), key=str)

    def _scored_steps(
        self, row: Row, remainder: PolyElement
    ) -> Iterator[tuple[int, Row]]:
        """Yields the rows with one more monomial and their remainders' distance.

        The determinant is linear in each column: adding a monomial m to the
        entry in column col changes it by Σᵢ shiftⁱ(m)⋅cofactorᵢ. The added
        monomial has the same degree as the entry.
        """
        for col in range(3):
            cofactors = self.cofactors(row, col)
            present = set(row[col].monoms())
            degree = sum(next(iter(present)))
            for exps, sign in itertools.product(self.monomials[degree], (1, -1)):
                if exps in present:
                    continue
                addend = self.monomial(exps, sign)
                change = sum(
                    (self.shift(addend, i) * c for i, c in enumerate(cofactors)),
                    self.ring.zero,
                )
                new_row = (*row[:col], row[col] + addend, *row[col + 1 :])
                yield _distance(remainder - change), new_row

    def run(self, rng: random.Random) -> Row | None:
        """Runs one randomized search. Returns None if it gets stuck."""
        if not self.starts:
            return None
        row = rng.choice(self.starts)
        while remainder := self.target - self.det(row):
            steps = _closest(self._scored_steps(row, remainder), _distance(remainder))
            if not steps:
                return None
            row = rng.choice(steps)
        return row


def poly_as_det3_matrix(
    poly: Expr,
    *,
    max_degree: int | None = None,
    max_runs: int = 20,
    seed: int = 0,
) -> tuple[Expr, Matrix]:
    """Finds a monomial k and a 3x3 matrix M such that poly = k⋅det(M).

    k is the gcd of the terms: the rational gcd of the coefficients times the
    common monomial factor, or 0 if poly is 0. The entries of M are sums of
    ±monomials of degree at most max_degree, its rows are cyclic shifts of each
    other, and its column degrees are non-increasing. Raises ValueError if poly
    is not a cyclic polynomial with rational coefficients, or no matrix is found
    in max_runs randomized runs.
    """
    poly = sympify(poly).expand()
    if poly == 0:
        return Rational(0), zeros(3, 3)
    if not poly.is_polynomial() or not all(
        c.is_Rational for c in poly.as_coefficients_dict().values()
    ):
        msg = f"{poly} is not a polynomial with rational coefficients"
        raise ValueError(msg)

    k, primitive = poly.as_content_primitive()
    index_shift = _index_shift(primitive.free_symbols)
    ring = PolyRing(sorted(index_shift, key=str), ZZ)
    target = ring.from_expr(primitive)
    common_exps = tuple(min(col) for col in zip(*target.monoms(), strict=True))
    monomial = ring.from_dict({common_exps: 1})
    k *= monomial.as_expr()
    target = target.exquo(monomial)
    if target.LC < 0:  # makes the search independent of the sign of poly
        k, target = -k, -target
    search = _Search(target, index_shift, max_degree)
    if search.shift(target) != target:
        msg = f"{poly} is not invariant under the index shift 1 → 2 → 3 → 1"
        raise ValueError(msg)

    rng = random.Random(seed)  # noqa: S311 (suspicious-non-cryptographic-random-usage)
    for _ in range(max_runs):
        row = search.run(rng)
        if row is not None:
            return k, Matrix([[e.as_expr() for e in r] for r in search.rows(row)])
    msg = f"cannot express {poly} as a 3x3 determinant"
    if constants := sorted(
        s.name for s, t in index_shift.items() if s == t and s.name[-1] in "123"
    ):
        msg += f" (treated as constants: {', '.join(constants)})"
    raise ValueError(msg)


def poly_as_det3(
    poly: Expr,
    *,
    max_degree: int | None = None,
    max_runs: int = 20,
    seed: int = 0,
) -> Expr:
    """Expresses poly as k⋅det(M) with an unevaluated 3x3 determinant.

    See poly_as_det3_matrix() for the details.
    """
    k, matrix = poly_as_det3_matrix(
        poly, max_degree=max_degree, max_runs=max_runs, seed=seed
    )
    return k * Determinant(matrix)


def normalize_det3(det: Expr) -> Expr:
    """Rewrites k⋅det(M) into a canonical form with the same value.

    k may be any expression. The steps:
    1. Negates rows and columns of M to minimize the number of minus signs,
       then the number of entries with a negative leading term.
    2. Sorts the columns by non-increasing degree, then the rows and columns by
       the variable indices of the entries in reading order (x1 before x2).
    3. Moves the sign change of the determinant into k.
    """
    if det == 0:
        return det
    k, determinant = det.as_independent(Determinant, as_Add=False)
    if not isinstance(determinant, Determinant):
        msg = f"{det} is not of the form k⋅det(M)"
        raise TypeError(msg)
    matrix = determinant.arg
    gens = sorted(matrix.free_symbols, key=str)
    if not gens:
        return det
    polys = [[Poly(e, *gens) for e in row] for row in matrix.tolist()]

    def minus_signs(signs: tuple[int, ...]) -> tuple[int, int]:
        """Minus signs in the matrix, then entries with a negative leading term.

        Row i is multiplied by signs[i], column j by signs[3 + j].
        """
        minus = [
            [c.is_negative != (signs[i] * signs[3 + j] < 0) for c in p.coeffs()]
            for i, row in enumerate(polys)
            for j, p in enumerate(row)
        ]
        return sum(map(sum, minus)), sum(m[0] for m in minus if m)

    # Negating all rows and columns is a no-op, so the last column stays as is.
    signs = min(
        ((*s, 1) for s in itertools.product((1, -1), repeat=5)), key=minus_signs
    )

    degrees = [max(polys[i][j].total_degree() for i in range(3)) for j in range(3)]
    keys = [[[tuple(-e for e in m) for m in p.monoms()] for p in row] for row in polys]

    def sorted_rows(cols: tuple[int, ...]) -> list[int]:
        return sorted(range(3), key=lambda i: [keys[i][j] for j in cols])

    def order_key(cols: tuple[int, ...]) -> tuple:
        rows = sorted_rows(cols)
        return [-degrees[j] for j in cols], [keys[i][j] for i in rows for j in cols]

    cols = min(itertools.permutations(range(3)), key=order_key)
    rows = sorted_rows(cols)

    signed = Matrix(3, 3, lambda i, j: signs[i] * signs[3 + j] * matrix[i, j])
    sign = math.prod(signs) * Permutation(rows).signature()
    sign *= Permutation(cols).signature()
    return k * sign * Determinant(signed.extract(rows, cols))
