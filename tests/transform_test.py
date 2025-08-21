from sympy import pi, simplify, symbols

from lib.central_conic import ConicCenter
from lib.circle import Circle
from lib.line import X_AXIS, HorizontalLine, LineBetween
from lib.matrix import ConicMatrix, IsNonZeroMultiple
from lib.transform import (
    Rotate,
    Scale,
    ScaleXY,
    TransformConic,
    TransformLine,
    Translate,
)

conic = ConicMatrix(*symbols("a,b,c,d,e,f"))


class TestTranslate:
    def test_translate_circle(self):
        x, y, r = symbols("x,y,r")
        dx, dy = symbols("dx,dy")
        transformation = Translate(dx, dy)
        circle = Circle((x, y), r)
        new_circle = TransformConic(circle, transformation)
        new_center_x, new_center_y = ConicCenter(new_circle)
        assert new_center_x == x + dx
        assert new_center_y == y + dy

    def test_translate_conic(self):
        center_x, center_y = ConicCenter(conic)
        dx, dy = symbols("dx,dy")
        transformation = Translate(dx, dy)
        new_conic = TransformConic(conic, transformation)
        new_center_x, new_center_y = ConicCenter(new_conic)
        assert dx == simplify(new_center_x - center_x)
        assert dy == simplify(new_center_y - center_y)


class TestRotate:
    def test_rotate_circle_around_center(self):
        x, y, r, theta = symbols("x,y,r,theta")
        circle = Circle((x, y), r)
        rotation = Rotate(theta, x, y)
        new_circle = TransformConic(circle, rotation)
        assert circle == simplify(new_circle)

    def test_rotation_around_point(self):
        x0, y0, theta = symbols("x0,y0,theta")
        rotation = Rotate(theta, x0, y0)
        rotation_sequence = Translate(x0, y0) * Rotate(theta) * Translate(-x0, -y0)
        assert simplify(rotation) == simplify(rotation_sequence)


class TestScale:
    def test_scaling_around_point(self):
        x0, y0, s = symbols("x0,y0,s")
        scaling = Scale(s, x0, y0)
        scaling_sequence = Translate(x0, y0) * Scale(s) * Translate(-x0, -y0)
        assert simplify(scaling) == simplify(scaling_sequence)

    def test_scaling_unevenly_around_point(self):
        x0, y0, sx, sy = symbols("x0,y0,sx,sy")
        scaling = ScaleXY(sx, sy, x0, y0)
        scaling_sequence = Translate(x0, y0) * ScaleXY(sx, sy) * Translate(-x0, -y0)
        assert simplify(scaling) == simplify(scaling_sequence)


class TestTransformLine:
    def test_rotate_line(self):
        rotated = TransformLine(X_AXIS, Rotate(pi / 4))
        assert IsNonZeroMultiple(rotated, LineBetween((0, 0), (1, 1)))

    def test_translate_line(self):
        translated = TransformLine(X_AXIS, Translate(1, 2))
        assert translated == HorizontalLine(2)
