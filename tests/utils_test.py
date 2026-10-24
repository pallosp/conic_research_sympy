from sympy import ImmutableMatrix, Matrix, Transpose

from tests.utils import are_projective_sets_equal


class TestAreProjectiveSetsEqual:
    def test_matrices(self):
        set1 = [Matrix([1, 0, 0]), Matrix([1, 2, 3])]
        set2 = [Matrix([2, 4, 6]), Matrix([-1, 0, 0])]
        assert are_projective_sets_equal(set1, set2)

    def test_lists_and_tuples(self):
        assert are_projective_sets_equal([[1, 0, 0], (1, 2, 3)], [(2, 4, 6), [3, 0, 0]])

    def test_matrix_expressions(self):
        column = Transpose(ImmutableMatrix([[1, 2, 3]]))
        assert are_projective_sets_equal(
            [Matrix([1, 0, 0]), column], [[2, 4, 6], [3, 0, 0]]
        )
        assert not are_projective_sets_equal([column], [[1, 2, 4]])

    def test_different_sizes(self):
        assert not are_projective_sets_equal([[1, 2, 3]], [[1, 2, 3], [1, 0, 0]])

    def test_multiplicity(self):
        assert not are_projective_sets_equal(
            [[1, 2, 3], [1, 2, 3]], [[1, 2, 3], [1, 0, 0]]
        )
