from sympy import Matrix
from sympy import symbols
from sympy import simplify

from lib.central_conic import ConicCenter
from lib.transform import TransformConic
from lib.transform import Translate

a, b, c, d, e, f = symbols("a,b,c,d,e,f")
conic = Matrix([[a, b, d], [b, c, e], [d, e, f]])


class TestTranslate:
    def test_translate_circle(self):
        center_x, center_y = ConicCenter(conic)
        dx, dy = symbols("dx,dy")
        transformation = Translate(dx, dy)
        new_conic = TransformConic(conic, transformation)
        new_center_x, new_center_y = ConicCenter(new_conic)
        assert dx == simplify(new_center_x - center_x)
        assert dy == simplify(new_center_y - center_y)
