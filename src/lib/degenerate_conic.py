from sympy import Matrix


def LinePair(line1: Matrix, line2: Matrix) -> Matrix:
    """Constructs a conic section from two projective lines."""
    assert line1.shape == (3, 1)
    assert line2.shape == (3, 1)
    return (line1 * line2.T + line2 * line1.T) / 2
