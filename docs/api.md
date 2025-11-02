# Table of Contents

* [matrix](#matrix)
  * [is\_nonzero\_multiple](#matrix.is_nonzero_multiple)
  * [is\_positive\_multiple](#matrix.is_positive_multiple)
  * [max\_eigenvalue](#matrix.max_eigenvalue)
  * [min\_eigenvalue](#matrix.min_eigenvalue)
  * [conic\_matrix](#matrix.conic_matrix)
  * [quadratic\_form](#matrix.quadratic_form)
  * [skew\_matrix](#matrix.skew_matrix)
  * [NonzeroCross](#matrix.NonzeroCross)
    * [eval](#matrix.NonzeroCross.eval)
  * [is\_real\_matrix](#matrix.is_real_matrix)
  * [is\_definite\_matrix](#matrix.is_definite_matrix)
  * [is\_full\_rank](#matrix.is_full_rank)
* [central\_conic](#central_conic)
  * [conic\_from\_foci\_and\_radius](#central_conic.conic_from_foci_and_radius)
  * [conic\_from\_center\_and\_points](#central_conic.conic_from_center_and_points)
  * [conic\_center](#central_conic.conic_center)
  * [semi\_axis\_lengths](#central_conic.semi_axis_lengths)
  * [primary\_radius](#central_conic.primary_radius)
  * [secondary\_radius](#central_conic.secondary_radius)
  * [linear\_eccentricity](#central_conic.linear_eccentricity)
  * [center\_to\_focus\_vector](#central_conic.center_to_focus_vector)
  * [center\_to\_vertex\_vector](#central_conic.center_to_vertex_vector)
  * [shrink\_conic\_to\_zero](#central_conic.shrink_conic_to_zero)
* [circle](#circle)
  * [circle](#circle.circle)
  * [circle\_radius](#circle.circle_radius)
  * [director\_circle](#circle.director_circle)
  * [UNIT\_CIRCLE](#circle.UNIT_CIRCLE)
  * [IMAGINARY\_UNIT\_CIRCLE](#circle.IMAGINARY_UNIT_CIRCLE)
* [line](#line)
  * [horizontal\_line](#line.horizontal_line)
  * [vertical\_line](#line.vertical_line)
  * [line\_between](#line.line_between)
  * [parallel\_line](#line.parallel_line)
  * [perpendicular\_line](#line.perpendicular_line)
  * [line\_through\_point](#line.line_through_point)
  * [angle\_bisector](#line.angle_bisector)
  * [perpendicular\_bisector](#line.perpendicular_bisector)
  * [line\_normal](#line.line_normal)
  * [are\_parallel](#line.are_parallel)
  * [are\_perpendicular](#line.are_perpendicular)
  * [IDEAL\_LINE](#line.IDEAL_LINE)
  * [X\_AXIS](#line.X_AXIS)
  * [Y\_AXIS](#line.Y_AXIS)
* [conic](#conic)
  * [conic\_from\_poly](#conic.conic_from_poly)
  * [conic\_through\_points](#conic.conic_through_points)
  * [conic\_from\_focus\_and\_directrix](#conic.conic_from_focus_and_directrix)
  * [eccentricity](#conic.eccentricity)
  * [focal\_axis\_direction](#conic.focal_axis_direction)
  * [focal\_axis](#conic.focal_axis)
  * [IdealPoints](#conic.IdealPoints)
    * [eval](#conic.IdealPoints.eval)
  * [projective\_conic\_center](#conic.projective_conic_center)
  * [pole\_point](#conic.pole_point)
  * [polar\_line](#conic.polar_line)
* [degenerate\_conic](#degenerate_conic)
  * [line\_pair\_conic](#degenerate_conic.line_pair_conic)
  * [double\_line\_conic](#degenerate_conic.double_line_conic)
  * [point\_conic](#degenerate_conic.point_conic)
  * [SplitToLines](#degenerate_conic.SplitToLines)
    * [eval](#degenerate_conic.SplitToLines.eval)
  * [ExtractPoint](#degenerate_conic.ExtractPoint)
    * [eval](#degenerate_conic.ExtractPoint.eval)
* [conic\_classes](#conic_classes)
  * [ConicNormFactor](#conic_classes.ConicNormFactor)
    * [eval](#conic_classes.ConicNormFactor.eval)
  * [is\_degenerate](#conic_classes.is_degenerate)
  * [is\_nondegenerate](#conic_classes.is_nondegenerate)
  * [is\_central\_conic](#conic_classes.is_central_conic)
  * [is\_finite\_conic](#conic_classes.is_finite_conic)
  * [is\_imaginary\_ellipse](#conic_classes.is_imaginary_ellipse)
  * [is\_ellipse](#conic_classes.is_ellipse)
  * [is\_circle](#conic_classes.is_circle)
  * [is\_parabola](#conic_classes.is_parabola)
  * [is\_hyperbola](#conic_classes.is_hyperbola)
  * [is\_circular](#conic_classes.is_circular)
  * [is\_line\_pair](#conic_classes.is_line_pair)
  * [is\_double\_line](#conic_classes.is_double_line)
  * [is\_point\_conic](#conic_classes.is_point_conic)
  * [is\_finite\_point\_conic](#conic_classes.is_finite_point_conic)
* [distance](#distance)
  * [point\_point\_distance](#distance.point_point_distance)
  * [point\_line\_distance](#distance.point_line_distance)
  * [parallel\_line\_distance](#distance.parallel_line_distance)
* [sympy\_utils](#sympy_utils)
  * [add\_eq](#sympy_utils.add_eq)
  * [sub\_eq](#sympy_utils.sub_eq)
  * [mul\_eq](#sympy_utils.mul_eq)
  * [div\_eq](#sympy_utils.div_eq)
  * [swap\_eq](#sympy_utils.swap_eq)
  * [eq\_chain](#sympy_utils.eq_chain)
* [point](#point)
  * [ORIGIN](#point.ORIGIN)
  * [ideal\_point](#point.ideal_point)
  * [ideal\_point\_on\_line](#point.ideal_point_on_line)
  * [point\_to\_xy](#point.point_to_xy)
  * [point\_to\_vec3](#point.point_to_vec3)
  * [centroid](#point.centroid)
  * [perpendicular\_foot](#point.perpendicular_foot)
* [transform](#transform)
  * [transform\_point](#transform.transform_point)
  * [transform\_line](#transform.transform_line)
  * [transform\_conic](#transform.transform_conic)
  * [translate](#transform.translate)
  * [rotate](#transform.rotate)
  * [reflect\_to\_line](#transform.reflect_to_line)
  * [scale\_xy](#transform.scale_xy)
  * [scale](#transform.scale)
  * [transformation\_from\_samples](#transform.transformation_from_samples)
* [parabola](#parabola)
  * [parabola\_directrix](#parabola.parabola_directrix)
  * [parabola\_focus](#parabola.parabola_focus)
  * [parabola\_vertex](#parabola.parabola_vertex)
  * [parabola\_direction](#parabola.parabola_direction)
  * [parabola\_axis](#parabola.parabola_axis)
  * [parabola\_focal\_parameter](#parabola.parabola_focal_parameter)
* [polar\_conic](#polar_conic)
  * [POLAR\_UNIT\_CIRCLE](#polar_conic.POLAR_UNIT_CIRCLE)
  * [point\_at\_angle](#polar_conic.point_at_angle)
  * [conic\_from\_polar\_matrix](#polar_conic.conic_from_polar_matrix)
* [incidence](#incidence)
  * [line\_contains\_point](#incidence.line_contains_point)
  * [conic\_contains\_point](#incidence.conic_contains_point)
  * [conic\_contains\_line](#incidence.conic_contains_line)
  * [are\_collinear](#incidence.are_collinear)
  * [are\_concurrent](#incidence.are_concurrent)
  * [are\_on\_same\_conic](#incidence.are_on_same_conic)
  * [are\_cocircular](#incidence.are_cocircular)
* [hyperbola](#hyperbola)
  * [hyperbola\_from\_foci\_and\_point](#hyperbola.hyperbola_from_foci_and_point)
* [ellipse](#ellipse)
  * [ellipse](#ellipse.ellipse)
  * [ellipse\_from\_foci\_and\_point](#ellipse.ellipse_from_foci_and_point)
  * [steiner\_ellipse](#ellipse.steiner_ellipse)
  * [steiner\_inellipse](#ellipse.steiner_inellipse)
* [intersection](#intersection)
  * [line\_x\_line](#intersection.line_x_line)
  * [conic\_x\_line](#intersection.conic_x_line)
* [transform\_classes](#transform_classes)
  * [is\_affine\_transform](#transform_classes.is_affine_transform)

<a id="matrix"></a>

# matrix

<a id="matrix.is_nonzero_multiple"></a>

#### is\_nonzero\_multiple

```python
def is_nonzero_multiple(m1: Matrix | Sequence[Expr],
                        m2: Matrix | Sequence[Expr]) -> bool | None
```

Tells whether two matrices are non-zero scalar multiples of each other.

Treats lists and tuples as column vectors. Returns `None` if undecidable.

<a id="matrix.is_positive_multiple"></a>

#### is\_positive\_multiple

```python
def is_positive_multiple(m1: Matrix | Sequence[Expr],
                         m2: Matrix | Sequence[Expr]) -> bool | None
```

Tells whether two matrices are positive scalar multiples of each other.

Treats lists and tuples as column vectors. Returns `None` if undecidable.

<a id="matrix.max_eigenvalue"></a>

#### max\_eigenvalue

```python
def max_eigenvalue(symmetric_matrix_2x2: Matrix) -> Expr
```

Returns the higher eigenvalue of a 2x2 symmetric matrix.

<a id="matrix.min_eigenvalue"></a>

#### min\_eigenvalue

```python
def min_eigenvalue(symmetric_matrix_2x2: Matrix) -> Expr
```

Returns the lower eigenvalue of a 2x2 symmetric matrix.

<a id="matrix.conic_matrix"></a>

#### conic\_matrix

```python
def conic_matrix(a: Expr, b: Expr, c: Expr, d: Expr, e: Expr,
                 f: Expr) -> Matrix
```

Builds a 3x3 symmetric conic matrix from its elements.

The conic equation looks like this:
```
        [a b d] [x]
[x y 1] [b c e] [y] = 0
        [d e f] [1]
```

or in expanded form `ax² + 2bxy + cy² + 2dx + 2ey + f = 0`

<a id="matrix.quadratic_form"></a>

#### quadratic\_form

```python
def quadratic_form(sym_matrix: Matrix, vector: Matrix) -> Expr
```

Computes the quadratic form for a n⨯n symmetric matrix and an n-element
column vector.

Use case: when `sym_matrix` and `vector` represent a conic and a projective
point, respectively, the quadratic form is zero iff the point is on the
conic.

*Formula*: `vᵀ·M·v`

<a id="matrix.skew_matrix"></a>

#### skew\_matrix

```python
def skew_matrix(vector3: Matrix) -> Matrix
```

Creates a skew-symmetric matrix from a 3D vector.

<a id="matrix.NonzeroCross"></a>

## NonzeroCross Objects

```python
class NonzeroCross(Function)
```

Finds a column and a row in a matrix whose intersection is a non-zero
element.

Returns an unevaluated `sympy.Function` if none of the elements can be
proven to be non-zero, or `nan` in case of a zero matrix.

<a id="matrix.NonzeroCross.eval"></a>

#### eval

```python
@classmethod
def eval(cls, matrix: Matrix) -> tuple[Matrix, Matrix] | NaN | None
```

Internal implementation. Call `NonZeroCross(matrix)` directly.

<a id="matrix.is_real_matrix"></a>

#### is\_real\_matrix

```python
def is_real_matrix(matrix: Matrix) -> bool | None
```

Checks if all elements of a matrix are real.

Returns a bool or `None` if undecidable.

<a id="matrix.is_definite_matrix"></a>

#### is\_definite\_matrix

```python
def is_definite_matrix(matrix: Matrix) -> bool | None
```

Checks if a real matrix is either positive or negative definite.

Returns a bool or `None` if the definitess can't be decided.

<a id="matrix.is_full_rank"></a>

#### is\_full\_rank

```python
def is_full_rank(matrix: Matrix,
                 *,
                 simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Tells whether a matrix has full rank.

Takes an optional `simplifier` callback that simplifies the determinant
before it gets compared to zero. Returns `None` if the result is
undecidable.

<a id="central_conic"></a>

# central\_conic

<a id="central_conic.conic_from_foci_and_radius"></a>

#### conic\_from\_foci\_and\_radius

```python
def conic_from_foci_and_radius(focus1: Matrix | Sequence[Expr],
                               focus2: Matrix | Sequence[Expr],
                               radius: Expr) -> Matrix
```

Computes the ellipse or hyperbola with the given focus points and
radius, i.e. center-vertex distance.

If `radius` is negative, takes its absolute value.

*Formula*:
[research/conic_from_foci_and_radius.py](../src/research/conic_from_foci_and_radius.py)

<a id="central_conic.conic_from_center_and_points"></a>

#### conic\_from\_center\_and\_points

```python
def conic_from_center_and_points(center: Matrix | Sequence[Expr],
                                 p1: Matrix | Sequence[Expr],
                                 p2: Matrix | Sequence[Expr],
                                 p3: Matrix | Sequence[Expr]) -> Matrix
```

Computes the conic section with the given center and perimeter points.

May return
 - an ellipse;
 - a hyperbola;
 - a parallel line pair;
 - zero matrix if the solution is ambiguous, which happens when some of the
   (`center`, `pᵢ`, `pⱼ`) triples are collinear.

*Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)

<a id="central_conic.conic_center"></a>

#### conic\_center

```python
def conic_center(conic: Matrix) -> Matrix
```

Computes the center point of a conic.

Returns the point's coordinates as a 2D column vector.

*Formula*: [research/conic_center.py](../src/research/conic_center.py)

<a id="central_conic.semi_axis_lengths"></a>

#### semi\_axis\_lengths

```python
def semi_axis_lengths(conic: Matrix) -> tuple[Expr, Expr]
```

Computes the semi-axis lengths of a conic in no specific order.

To get the semi-focal or semi-transverse axis length (semi-major /
semi-minor in case of ellipses), call
[principal_radius](#central_conic.principal_radius) or
[secondary_radius](#central_conic.secondary_radius), respectively.

*Formula*: [research/conic_radii.py](../src/research/conic_radii.py)

<a id="central_conic.primary_radius"></a>

#### primary\_radius

```python
def primary_radius(conic: Matrix) -> Expr
```

Computes the center-vertex distance of a conic.

This corresponds to the semi-major axis length of real ellipses. In case of
[imaginary ellipses](#conic_classification.is_imaginary_ellipse) however the
focal axis is the shorter one in terms of absolute value.

The returned value is:
 - a positive number for ellipses and hyperbolas;
 - infinity for parabolas;
 - an imaginary number for imaginary ellipses;
 - `nan` for ideal point conics;
 - 0 for the other degenerate conics.

<a id="central_conic.secondary_radius"></a>

#### secondary\_radius

```python
def secondary_radius(conic: Matrix) -> Expr
```

Computes the semi-conjugate axis length of a conic.

This corresponds to the semi-minor axis length of real ellipses. In case of
of imaginary ellipses however the conjugate axis is the longer one in terms
absolute value. Hyperbolas intersect their conjugate axis at complex points,
therefore the secondary radius will be a complex number.

The returned value is:
 - a positive number for ellipses;
 - infinity for parabolas;
 - an imaginary number for hyperbolas and imaginary ellipses;
 - `nan` for ideal point conics;
 - 0 for the other degenerate conics.

<a id="central_conic.linear_eccentricity"></a>

#### linear\_eccentricity

```python
def linear_eccentricity(conic: Matrix) -> Expr
```

Computes the linear eccentricity of a conic section.

The linear eccentricity is the distance between the center and a focus
point.

Special cases:
 - zero for circles and imaginary circles;
 - a real number for imaginary ellipses, since they still have real center
   and foci;
 - infinity for parabolas;
 - `nan` for parallel and coincident line pairs;
 - `nan` for conics containing the ideal line;
 - zero for all other degenerate conics.

*Formula*: `√|r₁²-r₂²|` where `r₁` and `r₂` denote the primary and secondary
radii of the conic (i.e., the semi-axis lengths in the case of an ellipse).

<a id="central_conic.center_to_focus_vector"></a>

#### center\_to\_focus\_vector

```python
def center_to_focus_vector(conic: Matrix) -> Matrix
```

Returns the 2D vector from a conic's center to one of its foci.

The opposite vector points to the other focus.

The function is only meaningful for
[central conics](#conic_classification.is_central_conic). The result vector
will contain infinite or `nan` elements when the conic lacks a finite
center.

<a id="central_conic.center_to_vertex_vector"></a>

#### center\_to\_vertex\_vector

```python
def center_to_vertex_vector(conic: Matrix) -> Matrix
```

Returns the 2D vector from a conic's center to one of its vertices.

The opposite vector points to the other vertex.

The function is only meaningful for
[central conics](#conic_classification.is_central_conic). The result vector
will contain infinite or `nan` elements when the conic lacks a finite
center.

<a id="central_conic.shrink_conic_to_zero"></a>

#### shrink\_conic\_to\_zero

```python
def shrink_conic_to_zero(conic: Matrix) -> Matrix
```

Scales a conic section from its center with a factor of zero.

Turns hyperbolas to line pair conics consisting of their asymptotes, and
ellipses to point conics.

This transformation is only meaningful for
[central conics](#conic_classification.is_central_conic): for other
conic types the result matrix will have infinite or `nan` elements.

*Formula*:
[research/scale_conic_from_center.py](../src/research/scale_conic_from_center.py

<a id="circle"></a>

# circle

<a id="circle.circle"></a>

#### circle

```python
def circle(center: Matrix | Sequence[Expr], radius: Expr) -> Matrix
```

Creates a circle from its center and radius.

<a id="circle.circle_radius"></a>

#### circle\_radius

```python
def circle_radius(circle: Matrix) -> Expr
```

Computes the radius of a circle conic.

The result is not specified if the conic matrix is not a circle.
The computation is based on
[research/director_circle.py](../src/research/director_circle.py).

<a id="circle.director_circle"></a>

#### director\_circle

```python
def director_circle(conic: Matrix) -> Matrix
```

Computes the director circle of a conic. It's also called orthoptic
circle or Fermat–Apollonius circle.

*Definition*: https://en.wikipedia.org/wiki/Director_circle<br>
*Formula*: [research/director_circle.py](../src/research/director_circle.py)

<a id="circle.UNIT_CIRCLE"></a>

#### UNIT\_CIRCLE

The circle at the origin with radius 1.

<a id="circle.IMAGINARY_UNIT_CIRCLE"></a>

#### IMAGINARY\_UNIT\_CIRCLE

The circle at the origin with radius 𝑖.

<a id="line"></a>

# line

<a id="line.horizontal_line"></a>

#### horizontal\_line

```python
def horizontal_line(y: Expr) -> Matrix
```

Constructs a horizontal line with the given y-coordinate.

<a id="line.vertical_line"></a>

#### vertical\_line

```python
def vertical_line(x: Expr) -> Matrix
```

Constructs a vertical line with the given x-coordinate.

<a id="line.line_between"></a>

#### line\_between

```python
def line_between(point1: Matrix | Sequence[Expr],
                 point2: Matrix | Sequence[Expr]) -> Matrix
```

Connects two projective points with a line.

<a id="line.parallel_line"></a>

#### parallel\_line

```python
def parallel_line(with_line: Matrix,
                  through_point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a line through a point parallel to a line.

<a id="line.perpendicular_line"></a>

#### perpendicular\_line

```python
def perpendicular_line(to_line: Matrix,
                       through_point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a line through a point perpendicular to a line.

<a id="line.line_through_point"></a>

#### line\_through\_point

```python
def line_through_point(point: Matrix | Sequence[Expr],
                       *,
                       direction: Matrix | Sequence[Expr] = None,
                       normal: Matrix | Sequence[Expr] = None) -> Matrix
```

Constructs a line through a point with the given direction.

The direction can be specified as
 - a 2D direction vector: `direction=(dx, dy)`
 - an ideal point on the line: `direction=(dx, dy, 0)`
 - a 2D normal vector: `normal=(nx, ny)`
 - an ideal point on the perpendicular line: `normal=(nx, ny, 0)`

<a id="line.angle_bisector"></a>

#### angle\_bisector

```python
def angle_bisector(line1: Matrix, line2: Matrix) -> Matrix
```

Constructs the angle bisector of two lines.

Chooses the bisector whose points substituted into the lines' equations have
the same sign. Negate one of the lines to get the other angle bisector.

Special cases:
 - AngleBisector(parallel finite lines, opposite direction) = center line
 - AngleBisector(coincident finite lines, same direction) = `[0, 0, 0]ᵀ`
 - AngleBisector(other parallel finite lines, same direction) = ideal line
 - AngleBisector(ideal line, ideal line) = `[0, 0, 0]ᵀ`
 - AngleBisector(ideal line, finite line) = ideal line

*Formula*: [research/angle_bisector.py](../src/research/angle_bisector.py)

<a id="line.perpendicular_bisector"></a>

#### perpendicular\_bisector

```python
def perpendicular_bisector(point1: Matrix | Sequence[Expr],
                           point2: Matrix | Sequence[Expr]) -> Matrix
```

Constructs the perpendicular bisector of two points.

<a id="line.line_normal"></a>

#### line\_normal

```python
def line_normal(line: Matrix,
                *,
                toward: Matrix | Sequence[Expr] = None) -> Matrix
```

Returns the normal vector of a line, represented as an ideal point.

By default, the normal is chosen such that `line·normal > 0`.
If a `toward` point is provided, the returned vector points toward it in the
projective sense: `line·normal` and `line·toward` have the same sign.
For an ordinary Euclidean point (given by its `(x, y)` coordinates), this
corresponds to the normal vector pointing toward the semiplane containing that
point.

Returns `[0, 0, 0]ᵀ` if the `toward` point lies on the line.

<a id="line.are_parallel"></a>

#### are\_parallel

```python
def are_parallel(line1: Matrix, line2: Matrix) -> bool | None
```

Tells whether line1 and line2 are parallel.

Returns `True` if they are parallel, `False` if not, `None` if undecidable.
Considers the ideal line parallel to everything.

<a id="line.are_perpendicular"></a>

#### are\_perpendicular

```python
def are_perpendicular(line1: Matrix, line2: Matrix) -> bool | None
```

Tells whether line1 and line2 are perpendicular.

Returns `True` if they are perpendicular, `False` if not, `None` if
undecidable. Considers the ideal line perpendicular to everything.

<a id="line.IDEAL_LINE"></a>

#### IDEAL\_LINE

The projective line at infinity.

<a id="line.X_AXIS"></a>

#### X\_AXIS

The horizontal line at y=0. Substituting a point at y>0 to the line's
equation will evaluate to a positive value.

<a id="line.Y_AXIS"></a>

#### Y\_AXIS

The vertical line at x=0. Substituting a point at x<0 to the line's equation
will evaluate to a positive value.

<a id="conic"></a>

# conic

<a id="conic.conic_from_poly"></a>

#### conic\_from\_poly

```python
def conic_from_poly(poly: Expr | Poly,
                    *,
                    x: Symbol = abc.x,
                    y: Symbol = abc.y) -> Matrix
```

Constructs a conic matrix from a two-variable quadratic polynomial.

Input: `ax² + bxy + cy² + dx + ey + f`

Output:
```
[a,   b/2, d/2]
[b/2, c,   e/2]
[d/2, e/2, f  ]
```

<a id="conic.conic_through_points"></a>

#### conic\_through\_points

```python
def conic_through_points(p1: Matrix | Sequence[Expr],
                         p2: Matrix | Sequence[Expr],
                         p3: Matrix | Sequence[Expr],
                         p4: Matrix | Sequence[Expr],
                         p5: Matrix | Sequence[Expr]) -> Matrix
```

Computes the conic that goes through the given points.

Returns the conic matrix, or a zero matrix if the result is ambiguous.
The result is unique if no 2 points coincide and no 4 points are collinear.
The result is non-degenerate if no 3 points are collinear.

*Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
section 10.1

<a id="conic.conic_from_focus_and_directrix"></a>

#### conic\_from\_focus\_and\_directrix

```python
def conic_from_focus_and_directrix(focus: Matrix | Sequence[Expr],
                                   directrix: Matrix,
                                   eccentricity: Expr) -> Matrix
```

Constructs a conic from its focus, directrix and eccentricity.

*Formula*:
[research/conic_from_focus_and_directrix.py](../src/research/conic_from_focus_and_directrix.py)

<a id="conic.eccentricity"></a>

#### eccentricity

```python
def eccentricity(conic: Matrix) -> Expr
```

Computes the eccentricity of a conic section.

The result is
 - 0 for circles and imaginary circles;
 - (0..1) for other real ellipses;
 - imaginary for other imaginary ellipses;
 - 1 for parabolas;
 - &gt;1 for hyperbolas, in particular √2 for rectangular hyperbolas.

In case of [non-degenerate](#conic_classification.is_nondegenerate)
[central conics](#conic_classification.is_central_conic), the eccentricity
equals to the ratio of the center-focus distance
([linear_eccentricity](#central_conic.linear_eccentricity))
and the center-vertex distance
([primary_radius](#central_conic.primary_radius)).

The eccentricity of finite point conics constructed by
[shrinking an ellipse to zero size](#central_conic.shrink_conic_to_zero)
equals to that of the original ellipse.

Crossing line pair conics have two different (generalized) focal axes, with
two different corresponding eccentricity values. Evaluate
`(eccentricity(conic), eccentricity(-conic))` to get both.

*Formula*:
https://en.wikipedia.org/wiki/Conic_section#Eccentricity_in_terms_of_coefficients<br>
*Own research*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

<a id="conic.focal_axis_direction"></a>

#### focal\_axis\_direction

```python
def focal_axis_direction(conic: Matrix) -> Matrix
```

Returns the ideal point representing the direction of a conic's focal axis.

Properties:
- The focal axis is treated as an undirected line; its angle to the
  horizontal lies in the (-π/2, π/2] interval. For the full direction of a
  parabola, use [parabola_direction](#parabola.parabola_direction) instead.
- Returns `[0, 0, 0]ᵀ` for circles and
  [circular conics](#conic_classification.is_circular).
- Point conics constructed by
  [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(ellipse)
  preserve the axis direction of the original real or imaginary ellipse.
- Line pair conics constructed by
  [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(hyperbola)
  have no such property.
- The focal axis of `line_pair(l1, l2)`, and `angle_bisector(l1, l2)` point
  to the same direction.

*Formula*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

<a id="conic.focal_axis"></a>

#### focal\_axis

```python
def focal_axis(conic: Matrix) -> Matrix
```

Returns the axis of symmetry going through conic's focus point(s).

Properties:
- Returns `[0, 0, 0]ᵀ` for circles and
  [circular conics](#conic_classification.is_circular).
- Point conics constructed by
  [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(ellipse)
  preserve the focal axis of the original real or imaginary ellipse.
- Line pair conics constructed by
  [shrink_conic_to_zero](#central_conic.shrink_conic_to_zero)(hyperbola)
  have no such property.
- The focal axis of `line_pair(l1, l2)`, and `angle_bisector(l1, l2)`
  coincide.

<a id="conic.IdealPoints"></a>

## IdealPoints Objects

```python
class IdealPoints(Function)
```

Computes the ideal points on a conic section.

Returns two points. Special cases:
 - For parabolas these are the same point.
 - For ellipses these are the complex conjugates of each other.
 - For symbolic conics returns an unevaluated `sympy.Function`.

<a id="conic.IdealPoints.eval"></a>

#### eval

```python
@classmethod
def eval(cls, conic: Matrix) -> tuple[Matrix, Matrix] | None
```

Internal implementation. Call `IdealPoints(conic)` directly.

<a id="conic.projective_conic_center"></a>

#### projective\_conic\_center

```python
def projective_conic_center(conic: Matrix) -> Matrix
```

Computes the generalized projective center of a conic.

It's equivalent to [conic_center](#central_conic.conic_center) (returns a
finite point) for
 - real and imaginary ellipses
 - hyperbolas
 - conics consisting of a single finite (Euclidean) point
 - crossing finite line pairs

For parabolas returns the ideal point on it
([proof](../src/research/parabola_center.py))

For other line pair conics and ideal point conics returns `[0, 0, 0]ᵀ`.

<a id="conic.pole_point"></a>

#### pole\_point

```python
def pole_point(conic: Matrix, polar_line: Matrix) -> Matrix
```

Computes the pole point of a conic with respect to the given polar line.

If the conic is degenerate, i.e. it factors into `l₁` and `l₂` real or
complex conjugate lines, the pole is
 - `[0, 0, 0]ᵀ` if `l₁`, `l₂`, and `polar_line` are concurrent;
 - the intersection of `l₁` and `l₂` otherwise.

*Pole / polar identity*: `conic * pole_point = polar_line`<br>
*Source*: https://en.wikipedia.org/wiki/Pole_and_polar#Calculating_the_pole_of_a_line

<a id="conic.polar_line"></a>

#### polar\_line

```python
def polar_line(conic: Matrix, pole_point: Matrix | Sequence[Expr]) -> Matrix
```

Computes the polar line of a conic with respect to the given pole point.

*Pole / polar identity*: `conic * pole_point = polar_line`
*Source*: https://en.wikipedia.org/wiki/Pole_and_polar#Calculating_the_pole_of_a_line

<a id="degenerate_conic"></a>

# degenerate\_conic

<a id="degenerate_conic.line_pair_conic"></a>

#### line\_pair\_conic

```python
def line_pair_conic(line1: Matrix, line2: Matrix) -> Matrix
```

Constructs a conic section from two projective lines.

<a id="degenerate_conic.double_line_conic"></a>

#### double\_line\_conic

```python
def double_line_conic(line: Matrix) -> Matrix
```

Constructs a degenerate conic consisting of two coincident lines.

<a id="degenerate_conic.point_conic"></a>

#### point\_conic

```python
def point_conic(point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a conic that degenerates to a single point.

Let the point's homogeneous coordinates be `(x, y, z)` and let
the variable point on the conic be `v = (X, Y, Z)`. The conic is
defined by the quadratic form `vᵀ C v = 0`.

For a conic consisting of just the given point, the condition
`x : y : z = X : Y : Z` must hold. Such conics can be expressed by the
equation

    λ₁(Y*z - Z*y)² + λ₂(Z*x - X*z)² + λ₃(X*y - Y*x)² = 0

where `λ₁`, `λ₂` and `λ₃` are either all positive or all negative.
This function returns the corresponding conic matrix for the choice
`λ₁ = λ₂ = λ₃ = -1`.

Hint: Use [ExtractPoint](#degenerate_conic.ExtractPoint) to recover the
point from the resulting conic.

<a id="degenerate_conic.SplitToLines"></a>

## SplitToLines Objects

```python
class SplitToLines(Function)
```

Splits a degenerate conic into two lines.

Special cases:
 - For non-degenerate conics the result is unspecified.
 - For point conics the lines will be complex conjugates.
 - For symbolic conics returns an unevaluated `sympy.Function`.

*Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
section 11.1

<a id="degenerate_conic.SplitToLines.eval"></a>

#### eval

```python
@classmethod
def eval(cls, conic: Matrix) -> tuple[Matrix, Matrix] | None
```

Internal implementation. Call `SplitToLines(conic)` directly.

<a id="degenerate_conic.ExtractPoint"></a>

## ExtractPoint Objects

```python
class ExtractPoint(Function)
```

Extracts the point from a point conic or the intersection of the
lines from a line pair conic.

Returns `[0, 0, 0]ᵀ` for double line pairs, or an unspecified 3D
column vector if the conic is not degenerate. May return an unevaluated
`sympy.Function` for symbolic conic matrices.

*Algorithm*:

Let the rows of the conic matrix be `r₁`, `r₂` and `r₃`. The point `p` is
on the conic iff `pᵀ C p = 0`, i.e. `p·(r₁·p, r₂·p, r₃·p) = 0`.

`p = r₁⨯r₂` is a solution, because `r₁·(r₁⨯r₂) = 0`, `r₂·(r₁⨯r₂) = 0` and
`r₃·(r₁⨯r₂) = det C = 0`. So are `p = r₂⨯r₃` and `p = r₃⨯r₁`. When the
conic matrix is a rank 2 matrix (point conic or non-coincident line pair),
at least one of these is a non-zero vector.

<a id="degenerate_conic.ExtractPoint.eval"></a>

#### eval

```python
@classmethod
def eval(cls, degenerate_conic: Matrix) -> Matrix | None
```

Internal implementation. Call `ExtractPoint(conic)` directly.

<a id="conic_classes"></a>

# conic\_classes

<a id="conic_classes.ConicNormFactor"></a>

## ConicNormFactor Objects

```python
class ConicNormFactor(Function)
```

Computes a normalization factor (±1) for a conic matrix `C`.

When `C` is multiplied by this factor, the resulting conic has the
following properties:

- *Non-degenerate conics*: the conic equation evaluates to a positive
  value at the focus point(s), i.e. `[fx fy 1]ᵀ C [fx fy 1] > 0`.
- *Point conics*: The conic equation evaluates to ≤0 for all finite
  points `[x, y, 1]ᵀ`.
- *Line-pair conics*: no preferred normalization exists; the factor is
  always 1.
- *Symbolic conics*: may return an unevaluated `sympy.Function` if the
  conic type or determinant sign cannot be determined.
- `conic.det() * ConicNormFactor(conic) == Abs(conic.det())` holds for
  all conic types.

<a id="conic_classes.ConicNormFactor.eval"></a>

#### eval

```python
@classmethod
def eval(cls, conic: Matrix) -> int | None
```

Internal implementation. Call `ConicNormFactor(conic)` directly.

<a id="conic_classes.is_degenerate"></a>

#### is\_degenerate

```python
def is_degenerate(
        conic: Matrix,
        *,
        simplifier: Callable[[Expr], Expr] = lambda expr: expr) -> bool | None
```

Tells whether the conic is degenerate.

Degenerate conics consist of a single projective point or a pair of
projective lines. The zero matrix is also considered degenerate.

Takes an optional `simplifier` callback that simplifies the conic matrix
determinant before it gets compared to zero. Returns `None` if the result
is undecidable.

<a id="conic_classes.is_nondegenerate"></a>

#### is\_nondegenerate

```python
def is_nondegenerate(
        conic: Matrix,
        *,
        simplifier: Callable[[Expr], Expr] = lambda expr: expr) -> bool | None
```

Tells whether the conic is non-degenerate.

Non-degenerate conics include real or imaginary ellipses, parabolas and
hyperbolas.

Takes an optional `simplifier` callback that simplifies the conic matrix
determinant before it gets compared to zero. Returns `None` if the result
is undecidable.

<a id="conic_classes.is_central_conic"></a>

#### is\_central\_conic

```python
def is_central_conic(
        conic: Matrix,
        *,
        simplifier: Callable[[Expr], Expr] = lambda expr: expr) -> bool | None
```

Tells whether a conic has a finite center of symmetry.

Takes an optional `simplifier` callback that simplifies the central
conicness polynomial before it gets compared to zero. Returns `None` if
the result is undecidable.

<a id="conic_classes.is_finite_conic"></a>

#### is\_finite\_conic

```python
def is_finite_conic(
        conic: Matrix,
        *,
        simplifier: Callable[[Expr], Expr] = factor) -> bool | None
```

Tells whether all points on the conic are finite.

Takes an optional `simplifier` callback that simplifies the finiteness
polynomial before it gets compared to zero. Returns `None` if the result
is undecidable.

<a id="conic_classes.is_imaginary_ellipse"></a>

#### is\_imaginary\_ellipse

```python
def is_imaginary_ellipse(conic: Matrix) -> bool | None
```

Tells whether the conic is an imaginary ellipse.

Imaginary ellipses have real center and focus points, but imaginary radii
and eccentricity. In addition, all solutions of their conic equations are
points with complex coordinates.

Returns `None` if undecidable.

<a id="conic_classes.is_ellipse"></a>

#### is\_ellipse

```python
def is_ellipse(conic: Matrix) -> bool | None
```

Tells whether the conic is an ellipse.

Returns `False` for
[imaginary ellipses](#conic_classification.is_imaginary_ellipse),
or `None` if the conic's type is undecidable.

<a id="conic_classes.is_circle"></a>

#### is\_circle

```python
def is_circle(conic: Matrix) -> bool | None
```

Tells whether the conic is a circle.

Returns `None` if undecidable.

<a id="conic_classes.is_parabola"></a>

#### is\_parabola

```python
def is_parabola(conic: Matrix) -> bool | None
```

Tells whether the conic is a parabola.

Returns `None` if undecidable.

<a id="conic_classes.is_hyperbola"></a>

#### is\_hyperbola

```python
def is_hyperbola(conic: Matrix) -> bool | None
```

Tells whether the conic is a hyperbola.

Returns `None` if undecidable.

<a id="conic_classes.is_circular"></a>

#### is\_circular

```python
def is_circular(conic: Matrix) -> bool | None
```

Tells whether there is a single center point around which the conic is
invariant under all rotations.

Circles, imaginary circles, zero-radius circles have such circular symmetry.
Double ideal lines are not considered circular. Returns `None` if undecidable.

<a id="conic_classes.is_line_pair"></a>

#### is\_line\_pair

```python
def is_line_pair(conic: Matrix) -> bool | None
```

Tells whether the conic is the union of two projective lines.

Returns `None` if undecidable.

<a id="conic_classes.is_double_line"></a>

#### is\_double\_line

```python
def is_double_line(conic: Matrix) -> bool | None
```

Tells whether the conic consists of two coincident projective lines.

Returns `None` if undecidable.

<a id="conic_classes.is_point_conic"></a>

#### is\_point\_conic

```python
def is_point_conic(conic: Matrix) -> bool | None
```

Tells whether the conic consists of a single projective point.

Returns `None` if undecidable.

A conic is a point conic iff it's degenerate and splits to two lines with
complex coordinates.

<a id="conic_classes.is_finite_point_conic"></a>

#### is\_finite\_point\_conic

```python
def is_finite_point_conic(conic: Matrix) -> bool | None
```

Tells whether the conic consists of a single finite (Euclidean) point.

Returns `None` if undecidable.

<a id="distance"></a>

# distance

<a id="distance.point_point_distance"></a>

#### point\_point\_distance

```python
def point_point_distance(point1: Matrix | Sequence[Expr],
                         point2: Matrix | Sequence[Expr]) -> Expr
```

Computes the signed distance between two points.

Special cases:
 - the distance between finite and ideal points is infinity
 - the distance between two ideal points is `nan`

<a id="distance.point_line_distance"></a>

#### point\_line\_distance

```python
def point_line_distance(point: Matrix | Sequence[Expr], line: Matrix) -> Expr
```

Computes the signed distance between a point and a line.

Special cases:
 - the distance between finite points and the ideal line is infinity
 - the distance between ideal points and any projective lines is `nan`

<a id="distance.parallel_line_distance"></a>

#### parallel\_line\_distance

```python
def parallel_line_distance(line1: Matrix, line2: Matrix) -> Expr
```

Computes the signed distance between two parallel lines.

The return value is positive if the two lines have the same direction,
negative otherwise.

Special cases:
- Returns infinity if one of the lines is the ideal line.
- Returns `nan` if both lines are ideal lines.
- Raises a `ValueError` if the lines are provably not parallel.
- Returns an unspecified value if the lines cross at a finite point, but
  Sympy cannot prove this fact.

<a id="sympy_utils"></a>

# sympy\_utils

<a id="sympy_utils.add_eq"></a>

#### add\_eq

```python
def add_eq(*eqs: Eq) -> Eq
```

Adds multiple sympy equations.

<a id="sympy_utils.sub_eq"></a>

#### sub\_eq

```python
def sub_eq(eq0: Eq, eq1: Expr | Eq) -> Eq
```

Subtracts a sympy equation or expression from another equation.

<a id="sympy_utils.mul_eq"></a>

#### mul\_eq

```python
def mul_eq(eq: Eq, factor: Expr | Eq) -> Eq
```

Multiplies a sympy equation by a factor or another equation.

<a id="sympy_utils.div_eq"></a>

#### div\_eq

```python
def div_eq(eq: Eq, denom: Expr) -> Eq
```

Divides a sympy equation by a denominator.

<a id="sympy_utils.swap_eq"></a>

#### swap\_eq

```python
def swap_eq(eq: Eq) -> Eq
```

Swaps the lhs and rhs of a sympy equation.

<a id="sympy_utils.eq_chain"></a>

#### eq\_chain

```python
def eq_chain(*expressions: Expr) -> Eq | Expr
```

Creates a chain of equations, i.e. `expr_1 = expr_2 = ... = expr_n`.

<a id="point"></a>

# point

<a id="point.ORIGIN"></a>

#### ORIGIN

The point at (0, 0)

<a id="point.ideal_point"></a>

#### ideal\_point

```python
def ideal_point(x: Expr, y: Expr) -> Matrix
```

Creates an ideal point at the given direction.

<a id="point.ideal_point_on_line"></a>

#### ideal\_point\_on\_line

```python
def ideal_point_on_line(line: Matrix) -> Matrix
```

Returns the coordinates of the ideal point on the line.

The first two coordinates specify the line's direction. The third one is
always zero.

If the line is the ideal line, returns `[0, 0, 0]ᵀ`.

<a id="point.point_to_xy"></a>

#### point\_to\_xy

```python
def point_to_xy(point: Matrix | Sequence[Expr]) -> tuple[Expr, Expr]
```

Computes the Euclidean coordinates of a projective point.

<a id="point.point_to_vec3"></a>

#### point\_to\_vec3

```python
def point_to_vec3(point: Matrix | Sequence[Expr]) -> Matrix
```

Computes the homogeneous coordinates of a projective point.

<a id="point.centroid"></a>

#### centroid

```python
def centroid(*points: Sequence[Expr]) -> Matrix
```

Computes the centroid of a set of points.

Returns the point's coordinates as a 2D column vector.

<a id="point.perpendicular_foot"></a>

#### perpendicular\_foot

```python
def perpendicular_foot(point: Matrix | Sequence[Expr], line: Matrix) -> Matrix
```

Computes the foot of the perpendicular through a point to a line.

Returns the foot point's coordinates as a 2D column vector.
Degenerates to `[nan, nan]ᵀ` when `point`, `line` or both are infinite.

*Formula*: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

<a id="transform"></a>

# transform

<a id="transform.transform_point"></a>

#### transform\_point

```python
def transform_point(point: Matrix | Sequence[Expr],
                    transformation: Matrix) -> Matrix
```

Applies a projective transformation to a projective point.

Returns a 3D column vector even if the input point is specified with only
two coordinates.

<a id="transform.transform_line"></a>

#### transform\_line

```python
def transform_line(line: Matrix, transformation: Matrix) -> Matrix
```

Applies a projective transformation to a projective line.

<a id="transform.transform_conic"></a>

#### transform\_conic

```python
def transform_conic(conic: Matrix, transformation: Matrix) -> Matrix
```

Applies a projective transformation to a conic.

<a id="transform.translate"></a>

#### translate

```python
def translate(dx: Expr, dy: Expr) -> Matrix
```

Computes the transformation matrix for a 2D translation.

<a id="transform.rotate"></a>

#### rotate

```python
def rotate(angle: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix
```

Computes the transformation matrix for a rotation around a point.

<a id="transform.reflect_to_line"></a>

#### reflect\_to\_line

```python
def reflect_to_line(axis: Matrix) -> Matrix
```

Computes the transformation matrix for a reflection to a line.

Returns a `nan` matrix if `axis` is the ideal line.

*Formula*:
[research/reflection_matrix.py](../src/research/reflection_matrix.py)

<a id="transform.scale_xy"></a>

#### scale\_xy

```python
def scale_xy(scale_x: Expr,
             scale_y: Expr,
             x0: Expr = 0,
             y0: Expr = 0) -> Matrix
```

Computes the projective transformation matrix for scaling along the x-
and y-axes.

<a id="transform.scale"></a>

#### scale

```python
def scale(scale: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix
```

Computes the projective transformation matrix for a uniform scaling
transformation.

<a id="transform.transformation_from_samples"></a>

#### transformation\_from\_samples

```python
def transformation_from_samples(
        source_points: Sequence[Matrix | Sequence[Expr]],
        target_points: Sequence[Matrix | Sequence[Expr]]) -> Matrix
```

Computes the transformation that maps one quadrilateral to another.

Takes 4 projective source and 4 projective target points. Both the source
and the target quadrilaterals must be non-degenerate. Returns a 3x3
projective transformation matrix.

*Research*:
[research/homography_from_samples.py](../src/research/homography_from_samples.py)

<a id="parabola"></a>

# parabola

<a id="parabola.parabola_directrix"></a>

#### parabola\_directrix

```python
def parabola_directrix(parabola: Matrix) -> Matrix
```

Computes the directrix of a parabola represented as a conic matrix.

Special cases for other conic types:
- Returns `[0, 0, 0]ᵀ` for coincident line pairs.
- Returns the ideal line for
  - non-coincident parallel line pairs;
  - conics consisting of one finite and one ideal line;
  - ideal point conics.
- Raises `ValueError` if the conic provably has 0 or 2 ideal points.
- Returns an unspecified 3D column vector in all other cases.

*Formula*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

<a id="parabola.parabola_focus"></a>

#### parabola\_focus

```python
def parabola_focus(parabola: Matrix) -> Matrix
```

Computes the focus of a parabola represented as a conic matrix.

Special cases for other conic types:
- Returns `[0, 0, 0]ᵀ` for
  - degenerate conics with one ideal point;
  - conics containing the ideal line.
- Raises `ValueError` if the conic provably has 0 or 2 ideal points.
- Returns an unspecified 3D column vector in all other cases.

*Formula*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

<a id="parabola.parabola_vertex"></a>

#### parabola\_vertex

```python
def parabola_vertex(parabola: Matrix) -> Matrix
```

Computes the parabola's vertex.

Returns the point's coordinates as a 2D column vector.

*Formula*: [research/parabola_vertex.py](../src/research/parabola_vertex.py)

<a id="parabola.parabola_direction"></a>

#### parabola\_direction

```python
def parabola_direction(parabola: Matrix) -> Matrix
```

Computes the direction of a parabola modulo 2π.

Unlike [focal_axis_direction](#conic.focal_axis_direction), which
determines the direction only modulo π, this function resolves the full
orientation.

Returns the ideal point representing the parabola’s direction. This point
coincides with both the ideal point lying on the parabola and its
[projective_conic_center](#conic.projective_conic_center).

<a id="parabola.parabola_axis"></a>

#### parabola\_axis

```python
def parabola_axis(parabola: Matrix) -> Matrix
```

Computes the parabola's focal axis line.

It's the polar line corresponding to the ideal point on the directrix.

Special cases for other conic types:
- Returns `[0, 0, 0]ᵀ` for
  - degenerate conics with one ideal point;
  - conics containing the ideal line.
- Raises `ValueError` if the conic provably has 0 or 2 ideal points.
- Returns an unspecified 3D column vector in all other cases.

<a id="parabola.parabola_focal_parameter"></a>

#### parabola\_focal\_parameter

```python
def parabola_focal_parameter(parabola: Matrix) -> Expr
```

Computes the parabola's focus-directrix distance.

*Formula*: [research/focal_parameter.py](../src/research/focal_parameter.py)

<a id="polar_conic"></a>

# polar\_conic

Utilities for conic sections in polar parametric form.

Polar conic sections are represented as parametric curves with angle θ, where
each point is the image of a point on the unit circle (`x² + y² = 1`) under a
projective transformation:

```
       [a b c]   [cos θ]
C(θ) = [d e f] * [sin θ]
       [g h i]   [  1  ]
```

<a id="polar_conic.POLAR_UNIT_CIRCLE"></a>

#### POLAR\_UNIT\_CIRCLE

The circle at the origin with radius 1, in polar matrix form.

<a id="polar_conic.point_at_angle"></a>

#### point\_at\_angle

```python
def point_at_angle(polar: Matrix, theta: Expr) -> Matrix
```

Computes the coordinates of the projective point on a polar conic
corresponding to a certain angle.

<a id="polar_conic.conic_from_polar_matrix"></a>

#### conic\_from\_polar\_matrix

```python
def conic_from_polar_matrix(polar: Matrix) -> Matrix
```

Transforms a conic from polar to quadratic form.

The algorithm is essentially applying the polar matrix as a projective
transformation on the unit circle.

<a id="incidence"></a>

# incidence

<a id="incidence.line_contains_point"></a>

#### line\_contains\_point

```python
def line_contains_point(
        line: Matrix,
        point: Matrix | Sequence[Expr],
        *,
        simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Tells whether `point` is on `line`.

Takes an optional `simplifier` callback that simplifies the incidence
polynomial before it gets compared to zero. Returns `None` if the result is
undecidable.

<a id="incidence.conic_contains_point"></a>

#### conic\_contains\_point

```python
def conic_contains_point(
        conic: Matrix,
        point: Matrix | Sequence[Expr],
        *,
        simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Checks if a point lies on a conic.

Takes an optional `simplifier` callback that simplifies the incidence
polynomial before it gets compared to zero. Returns `None` if the result is
undecidable.

<a id="incidence.conic_contains_line"></a>

#### conic\_contains\_line

```python
def conic_contains_line(
        conic: Matrix,
        line: Matrix,
        *,
        simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Checks if a line lies on a conic.

Takes an optional `simplifier` callback that simplifies the elements of the
containment matrix before it gets compared to zero. Returns `None` if the
result is undecidable.

*Formula*:
[research/conic_line_containment.py](../src/research/conic_line_containment.py)

<a id="incidence.are_collinear"></a>

#### are\_collinear

```python
def are_collinear(points: Sequence[Matrix],
                  *,
                  simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Tells whether n points are collinear.

Takes an optional `simplifier` callback that simplifies the collinearity
polynomial before it gets compared to zero. Returns `None` if undecidable.

<a id="incidence.are_concurrent"></a>

#### are\_concurrent

```python
def are_concurrent(lines: Sequence[Matrix],
                   *,
                   simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Tells whether n lines are concurrent, i.e. go through the same point.

Takes an optional callback that simplifies the concurrence polynomial
before it gets compared to zero. Returns `None` if undecidable.

Leverages the projective point-line duality, and uses the collinearity
formula described at
https://en.wikipedia.org/wiki/Incidence_(geometry)#Collinearity

<a id="incidence.are_on_same_conic"></a>

#### are\_on\_same\_conic

```python
def are_on_same_conic(
        points: Sequence[Matrix],
        *,
        simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Tells whether up to 6 points lie on the same conic section.

Takes an optional `simplifier` callback that simplifies the incidence
polynomial before it gets compared to zero. Returns `None` if undecidable.

*Formula*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
section 10.2 (Conics and Cross-Ratios)

<a id="incidence.are_cocircular"></a>

#### are\_cocircular

```python
def are_cocircular(points: Sequence[Matrix],
                   *,
                   simplifier: Callable[[Expr], Expr] = expand) -> bool | None
```

Tells whether n points lie on the same circle.

Note that four collinear points, as well as three collinear points and an
arbitrary ideal point are also considered cocircular.

Takes an optional `simplifier` callback that simplifies the cocircularity
determinant before it gets compared to zero. Returns `None` if undecidable.

*Formula*: [Measuring cocircularity in a point set](
https://egc23.web.uah.es/wp-content/uploads/2023/06/EGC23_paper_20.pdf)<br>
*Own research*:
[research/cocircularity.py](../src/research/cocircularity.py)

<a id="hyperbola"></a>

# hyperbola

<a id="hyperbola.hyperbola_from_foci_and_point"></a>

#### hyperbola\_from\_foci\_and\_point

```python
def hyperbola_from_foci_and_point(focus1: Matrix | Sequence[Expr],
                                  focus2: Matrix | Sequence[Expr],
                                  point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a hyperbola from its focus points and an incident point.

*Formula*:
[research/conic_from_foci_and_radius.py](../src/research/conic_from_foci_and_radius.py)

<a id="ellipse"></a>

# ellipse

<a id="ellipse.ellipse"></a>

#### ellipse

```python
def ellipse(center: Matrix | Sequence[Expr],
            r1: Expr,
            r2: Expr,
            *,
            r1_angle: Expr = None,
            r1_direction: Expr = None) -> Matrix
```

Constructs an ellipse from its center, radii, and the either the
direction vector of the first radius or its angle to horizontal.

*Formula*:
[research/ellipse_from_params.py](../src/research/ellipse_from_params.py)

<a id="ellipse.ellipse_from_foci_and_point"></a>

#### ellipse\_from\_foci\_and\_point

```python
def ellipse_from_foci_and_point(focus1: Matrix | Sequence[Expr],
                                focus2: Matrix | Sequence[Expr],
                                point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs an ellipse from its focus points and an incident point.

Special cases:
 - If all three points conicide, returns a zero matrix.
 - If the foci coincide but the point doesn't the result is a circle.
 - If the points are collinear but the foci don't coincide, returns a
   coincident line conic.
 - If any of the points are ideal points, returns a matrix that contains
   `nan` elements.

*Formula*:
[research/ellipse_from_foci_and_point.py](../src/research/ellipse_from_foci_and_point.py)

<a id="ellipse.steiner_ellipse"></a>

#### steiner\_ellipse

```python
def steiner_ellipse(point1: Matrix | Sequence[Expr],
                    point2: Matrix | Sequence[Expr],
                    point3: Matrix | Sequence[Expr]) -> Matrix
```

Constructs the Steiner circumellipse for the given points.

The ellipse goes through the three points and is centered at the triangle's
centroid.

*Definition*: https://en.wikipedia.org/wiki/Steiner_ellipse<br>
*Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)

<a id="ellipse.steiner_inellipse"></a>

#### steiner\_inellipse

```python
def steiner_inellipse(point1: Matrix | Sequence[Expr],
                      point2: Matrix | Sequence[Expr],
                      point3: Matrix | Sequence[Expr]) -> Matrix
```

Computes the Steiner inellipse for the given points.

The ellipse is centered at the triangle's centroid, and is tangent to the
triangle's sides at their midpoints.

*Definition*: https://en.wikipedia.org/wiki/Steiner_inellipse<br>
*Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)

<a id="intersection"></a>

# intersection

<a id="intersection.line_x_line"></a>

#### line\_x\_line

```python
def line_x_line(line1: Matrix, line2: Matrix) -> Matrix
```

Computes the intersection of two lines.

Returns an ideal point if the lines are parallel, or `[0, 0, 0]ᵀ` if they
coincide.

<a id="intersection.conic_x_line"></a>

#### conic\_x\_line

```python
def conic_x_line(
    conic: Matrix, line: Matrix
) -> tuple[Matrix | Sequence[Expr], Matrix | Sequence[Expr]] | NaN
```

Intersects a conic with a line. Returns two points.

Special cases:
 - The intersection points coincide if the line is tangent to the conic.
 - They are complex conjugates if the line doesn't intersect the conic at
   a real point.
 - Returns an unevaluated `sympy.Function` for symbolic conics.
 - Returns `None` if the conic contains the entire line.

*Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
section 11.3

<a id="transform_classes"></a>

# transform\_classes

<a id="transform_classes.is_affine_transform"></a>

#### is\_affine\_transform

```python
def is_affine_transform(transformation: Matrix) -> bool | None
```

Tells whether a transformation matrix is an affine transformation.

Returns None if undecidable.

