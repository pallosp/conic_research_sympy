from sympy import Matrix, expand, nan, pi, simplify, symbols

from lib.central_conic import conic_center
from lib.circle import circle
from lib.line import IDEAL_LINE, X_AXIS, horizontal_line, line_between
from lib.matrix import conic_matrix, is_nonzero_multiple
from lib.polar_conic import conic_from_polar_matrix, point_at_angle
from lib.transform import (
    homography_from_samples,
    reflect_to_line,
    rotate,
    scale,
    scale_xy,
    transform_conic,
    transform_line,
    transform_point,
    transform_polar_conic,
    translate,
)

conic = conic_matrix(*symbols("a,b,c,d,e,f"))


class TestTranslate:
    def test_translate_circle(self):
        x, y, r = symbols("x,y,r")
        dx, dy = symbols("dx,dy")
        transformation = translate((dx, dy))
        orig_circle = circle((x, y), r)
        translated_circle = transform_conic(orig_circle, transformation)
        new_center_x, new_center_y = conic_center(translated_circle)
        assert new_center_x == x + dx
        assert new_center_y == y + dy

    def test_translate_conic(self):
        center_x, center_y = conic_center(conic)
        dx, dy = symbols("dx,dy")
        transformation = translate((dx, dy))
        new_conic = transform_conic(conic, transformation)
        new_center_x, new_center_y = conic_center(new_conic)
        assert dx == simplify(new_center_x - center_x)
        assert dy == simplify(new_center_y - center_y)


class TestRotate:
    def test_rotate_circle_around_center(self):
        center = symbols("x y")
        r, angle = symbols("r theta")
        orig_circle = circle(center, r)
        rotation = rotate(angle, center)
        rotated_circle = transform_conic(orig_circle, rotation)
        assert orig_circle == simplify(rotated_circle)

    def test_rotation_around_point(self):
        angle = symbols("theta")
        center = Matrix(symbols("x y"))
        rotation = rotate(angle, center)
        rotation_sequence = translate(center) * rotate(angle) * translate(-center)
        assert simplify(rotation) == simplify(rotation_sequence)


class TestReflection:
    def test_reflect_to_finite_line(self):
        axis = Matrix([1, 1, -4])  # x+y=4
        point = Matrix([3, 0, 1])
        assert reflect_to_line(axis) * point == Matrix([4, 1, 1])

    def test_reflect_to_ideal_line(self):
        reflection = reflect_to_line(IDEAL_LINE)
        assert reflection == Matrix.ones(3, 3) * nan


class TestScale:
    def test_scaling_around_point(self):
        x0, y0, s = symbols("x0,y0,s")
        scaling = scale(s, center=(x0, y0))
        scaling_sequence = translate((x0, y0)) * scale(s) * translate((-x0, -y0))
        assert simplify(scaling) == simplify(scaling_sequence)

    def test_scaling_unevenly_around_point(self):
        sx, sy = symbols("sx sy")
        center = Matrix(symbols("x0 y0"))
        scaling = scale_xy(sx, sy, center)
        scaling_sequence = translate(center) * scale_xy(sx, sy) * translate(-center)
        assert simplify(scaling) == simplify(scaling_sequence)


class TestTransformPoint:
    def test_translate(self):
        assert transform_point((1, 2), translate((3, 5))) == Matrix([4, 7, 1])
        assert transform_point((1, 2, 0), translate((3, 5))) == Matrix([1, 2, 0])


class TestTransformLine:
    def test_rotate_line(self):
        rotated = transform_line(X_AXIS, rotate(pi / 4))
        assert is_nonzero_multiple(rotated, line_between((0, 0), (1, 1)))

    def test_translate_line(self):
        translated = transform_line(X_AXIS, translate((1, 2)))
        assert translated == horizontal_line(2)


class TestTransformPolarConic:
    def test_symbolic_transformation(self):
        conic = Matrix(3, 3, symbols("a b c d e f g h i"))
        transformation = Matrix(3, 3, symbols("A B C D E F G H I"))
        transformed_conic = transform_polar_conic(conic, transformation)

        phi = symbols("phi")
        point = point_at_angle(conic, phi)
        transformed_point = transform_point(point, transformation)
        point_on_transformed_conic = point_at_angle(transformed_conic, phi)

        assert expand(transformed_point) == expand(point_on_transformed_conic)

    def test_scaled_polar_conic_shape(self):
        polar_conic = Matrix(3, 3, symbols("a b c d e f g h i"))
        cartesian = conic_from_polar_matrix(polar_conic)

        scaling = scale_xy(2, 3, center=(4, 5))
        scaled_polar = transform_polar_conic(polar_conic, scaling)
        scaled_cartesian = transform_conic(cartesian, scaling)

        assert expand(conic_from_polar_matrix(scaled_polar)) == expand(scaled_cartesian)


class TestHomographyFromSamples:
    def test_rotate_90(self):
        source = ((1, 0), (0, 1), (-1, 0), (0, -1))
        target = ((0, 1), (-1, 0), (0, -1), (1, 0))
        transform = homography_from_samples(source, target)
        assert transform == rotate(pi / 2)

    def test_translate(self):
        source = ((0, 0), (1, 0), (0, 1), (1, 1))
        target = ((2, 1), (3, 1), (2, 2), (3, 2))
        transform = homography_from_samples(source, target)
        assert transform == translate((2, 1))

    def test_translate_homogeneous_coordinates(self):
        source = ((0, 0), (-2, 0, -2), (0, 1), (1, 1))
        target = ((-4, -2, -2), (3, 1), (2, 2), (3, 2))
        transform = homography_from_samples(source, target)
        assert transform == translate((2, 1))

    def test_circle_to_hyperbola(self):
        circle = ((1, 0), (0, 1), (-1, 0), (0, -1))
        hyperbola = ((1, 0, 1), (1, 1, 0), (-1, 0, 1), (1, -1, 0))
        transform = homography_from_samples(circle, hyperbola)
        for source, expected in zip(circle, hyperbola, strict=True):
            transformed = transform_point(source, transform)
            assert is_nonzero_multiple(transformed, expected)
        assert transform == Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
