# Table of Contents

* [matrix](#matrix)
  * [IsNonZeroMultiple](#matrix.IsNonZeroMultiple)
  * [MaxEigenvalue](#matrix.MaxEigenvalue)
  * [MinEigenvalue](#matrix.MinEigenvalue)
  * [ConicMatrix](#matrix.ConicMatrix)
  * [QuadraticForm](#matrix.QuadraticForm)
  * [SkewMatrix](#matrix.SkewMatrix)
  * [NonZeroCross](#matrix.NonZeroCross)
    * [eval](#matrix.NonZeroCross.eval)
  * [IsRealMatrix](#matrix.IsRealMatrix)
  * [IsDefinite](#matrix.IsDefinite)
* [central\_conic](#central_conic)
  * [ConicFromCenterAndPoints](#central_conic.ConicFromCenterAndPoints)
  * [ConicCenter](#central_conic.ConicCenter)
  * [SemiAxisLengths](#central_conic.SemiAxisLengths)
  * [SemiMajorAxis](#central_conic.SemiMajorAxis)
    * [eval](#central_conic.SemiMajorAxis.eval)
  * [SemiMinorAxis](#central_conic.SemiMinorAxis)
    * [eval](#central_conic.SemiMinorAxis.eval)
* [circle](#circle)
  * [Circle](#circle.Circle)
  * [CircleRadius](#circle.CircleRadius)
  * [DirectorCircle](#circle.DirectorCircle)
  * [UNIT\_CIRCLE](#circle.UNIT_CIRCLE)
  * [COMPLEX\_UNIT\_CIRCLE](#circle.COMPLEX_UNIT_CIRCLE)
* [line](#line)
  * [LineContainsPoint](#line.LineContainsPoint)
  * [HorizontalLine](#line.HorizontalLine)
  * [VerticalLine](#line.VerticalLine)
  * [LineBetween](#line.LineBetween)
  * [ParallelLine](#line.ParallelLine)
  * [PerpendicularLine](#line.PerpendicularLine)
  * [LineThroughPoint](#line.LineThroughPoint)
  * [AngleBisector](#line.AngleBisector)
  * [PerpendicularBisector](#line.PerpendicularBisector)
  * [AreParallel](#line.AreParallel)
  * [ArePerpendicular](#line.ArePerpendicular)
  * [IDEAL\_LINE](#line.IDEAL_LINE)
  * [X\_AXIS](#line.X_AXIS)
  * [Y\_AXIS](#line.Y_AXIS)
* [conic](#conic)
  * [ConicFromPoly](#conic.ConicFromPoly)
  * [ConicThroughPoints](#conic.ConicThroughPoints)
  * [ConicFromFocusAndDirectrix](#conic.ConicFromFocusAndDirectrix)
  * [Eccentricity](#conic.Eccentricity)
  * [AxisDirection](#conic.AxisDirection)
  * [IdealPoints](#conic.IdealPoints)
    * [eval](#conic.IdealPoints.eval)
  * [ProjectiveConicCenter](#conic.ProjectiveConicCenter)
  * [PolePoint](#conic.PolePoint)
  * [PolarLine](#conic.PolarLine)
  * [ConicContainsPoint](#conic.ConicContainsPoint)
  * [ConicContainsLine](#conic.ConicContainsLine)
* [degenerate\_conic](#degenerate_conic)
  * [LinePair](#degenerate_conic.LinePair)
  * [PointConic](#degenerate_conic.PointConic)
  * [SplitToLines](#degenerate_conic.SplitToLines)
    * [eval](#degenerate_conic.SplitToLines.eval)
  * [ExtractPoint](#degenerate_conic.ExtractPoint)
    * [eval](#degenerate_conic.ExtractPoint.eval)
* [distance](#distance)
  * [PointPointDistance](#distance.PointPointDistance)
  * [PointLineDistance](#distance.PointLineDistance)
* [sympy\_utils](#sympy_utils)
  * [AddEq](#sympy_utils.AddEq)
  * [SubEq](#sympy_utils.SubEq)
  * [MulEq](#sympy_utils.MulEq)
  * [DivEq](#sympy_utils.DivEq)
  * [SwapEq](#sympy_utils.SwapEq)
  * [EqChain](#sympy_utils.EqChain)
  * [FactorRadicals](#sympy_utils.FactorRadicals)
* [conic\_classification](#conic_classification)
  * [IsDegenerate](#conic_classification.IsDegenerate)
  * [IsNonDegenerate](#conic_classification.IsNonDegenerate)
  * [IsFiniteConic](#conic_classification.IsFiniteConic)
  * [IsComplexEllipse](#conic_classification.IsComplexEllipse)
  * [IsEllipse](#conic_classification.IsEllipse)
  * [IsCircle](#conic_classification.IsCircle)
  * [IsParabola](#conic_classification.IsParabola)
  * [IsHyperbola](#conic_classification.IsHyperbola)
  * [IsCircular](#conic_classification.IsCircular)
  * [IsLinePair](#conic_classification.IsLinePair)
  * [IsDoubleLine](#conic_classification.IsDoubleLine)
  * [IsPointConic](#conic_classification.IsPointConic)
  * [IsFinitePointConic](#conic_classification.IsFinitePointConic)
* [point](#point)
  * [ORIGIN](#point.ORIGIN)
  * [IdealPoint](#point.IdealPoint)
  * [IdealPointOnLine](#point.IdealPointOnLine)
  * [PointToXY](#point.PointToXY)
  * [PointToVec3](#point.PointToVec3)
  * [Centroid](#point.Centroid)
  * [PerpendicularFoot](#point.PerpendicularFoot)
* [transform](#transform)
  * [TransformConic](#transform.TransformConic)
  * [TransformLine](#transform.TransformLine)
  * [Translate](#transform.Translate)
  * [Rotate](#transform.Rotate)
  * [ScaleXY](#transform.ScaleXY)
  * [Scale](#transform.Scale)
* [parabola](#parabola)
  * [ParabolaDirectrix](#parabola.ParabolaDirectrix)
  * [ParabolaFocus](#parabola.ParabolaFocus)
* [polar\_conic](#polar_conic)
  * [POLAR\_UNIT\_CIRCLE](#polar_conic.POLAR_UNIT_CIRCLE)
  * [PointAtAngle](#polar_conic.PointAtAngle)
  * [ConicFromPolarMatrix](#polar_conic.ConicFromPolarMatrix)
* [ellipse](#ellipse)
  * [Ellipse](#ellipse.Ellipse)
  * [EllipseFromFociAndPoint](#ellipse.EllipseFromFociAndPoint)
  * [SteinerEllipse](#ellipse.SteinerEllipse)
  * [SteinerInellipse](#ellipse.SteinerInellipse)
* [intersection](#intersection)
  * [LineXLine](#intersection.LineXLine)
  * [ConicXLine](#intersection.ConicXLine)

<a id="matrix"></a>

# matrix

<a id="matrix.IsNonZeroMultiple"></a>

#### IsNonZeroMultiple

```python
def IsNonZeroMultiple(m1: Matrix | Sequence[Expr],
                      m2: Matrix | Sequence[Expr]) -> bool | None
```

Tells whether two matrices are non-zero scalar multiples of each other.

Treats lists and tuples as column vectors. Returns None if undecidable.

<a id="matrix.MaxEigenvalue"></a>

#### MaxEigenvalue

```python
def MaxEigenvalue(symmetric_matrix_2x2: Matrix) -> Expr
```

Returns the higher eigenvalue of a 2x2 symmetric matrix.

<a id="matrix.MinEigenvalue"></a>

#### MinEigenvalue

```python
def MinEigenvalue(symmetric_matrix_2x2: Matrix) -> Expr
```

Returns the lower eigenvalue of a 2x2 symmetric matrix.

<a id="matrix.ConicMatrix"></a>

#### ConicMatrix

```python
def ConicMatrix(a: Expr, b: Expr, c: Expr, d: Expr, e: Expr,
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

<a id="matrix.QuadraticForm"></a>

#### QuadraticForm

```python
def QuadraticForm(sym_matrix: Matrix, vector: Matrix) -> Expr
```

Computes the quadratic form for a nxn symmetric matrix and an n-element
column vector.

Use case: when `sym_matrix` and `vector` represent a conic and a projective
point, respectively, the quadratic form is zero iff the point is on the
conic.

*Formula*: `vᵀ·M·v`

<a id="matrix.SkewMatrix"></a>

#### SkewMatrix

```python
def SkewMatrix(vector3: Matrix) -> Matrix
```

Creates a skew-symmetric matrix from a 3d vector.

<a id="matrix.NonZeroCross"></a>

## NonZeroCross Objects

```python
class NonZeroCross(Function)
```

Finds a column and a row in a matrix whose intersection is a non-zero
element.

Returns an unevaluated `sympy.Function` if none of the elements can be
proven to be non-zero, or `nan` in case of a zero matrix.

<a id="matrix.NonZeroCross.eval"></a>

#### eval

```python
@classmethod
def eval(cls, matrix: Matrix) -> tuple[Matrix, Matrix] | NaN | None
```

Internal implementation. Call `NonZeroCross(matrix)` directly.

<a id="matrix.IsRealMatrix"></a>

#### IsRealMatrix

```python
def IsRealMatrix(matrix: Matrix) -> bool | None
```

Checks if all elements of a matrix are real.

Returns a bool or `None` if undecidable.

<a id="matrix.IsDefinite"></a>

#### IsDefinite

```python
def IsDefinite(matrix: Matrix) -> bool | None
```

Checks if a real matrix is either positive or negative definite.

Returns a bool or `None` if the definitess can't be decided.

<a id="central_conic"></a>

# central\_conic

<a id="central_conic.ConicFromCenterAndPoints"></a>

#### ConicFromCenterAndPoints

```python
def ConicFromCenterAndPoints(center: Matrix | Sequence[Expr],
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

<a id="central_conic.ConicCenter"></a>

#### ConicCenter

```python
def ConicCenter(conic: Matrix) -> tuple[Expr, Expr]
```

Computes the center point of a conic.

*Formula*: [research/conic_center.py](../src/research/conic_center.py)

<a id="central_conic.SemiAxisLengths"></a>

#### SemiAxisLengths

```python
def SemiAxisLengths(conic: Matrix) -> tuple[Expr, Expr]
```

Computes the semi-axis lengths of a conic.

*Formula*: [research/conic_radii.py](../src/research/conic_radii.py)

<a id="central_conic.SemiMajorAxis"></a>

## SemiMajorAxis Objects

```python
class SemiMajorAxis(Function)
```

Computes the semi-major axis length i.e. the center-vertex distance of
a conic.

The returned value is:
 - a positive number for ellipses and hyperbolas;
 - infinity for parabolas;
 - an imaginary number for complex ellipses;
 - nan for ideal point conics;
 - zero for the other degenerate conics.

Returns an unevaluated `sympy.Function` if we can't tell which axis is
longer.

<a id="central_conic.SemiMajorAxis.eval"></a>

#### eval

```python
@classmethod
def eval(cls, conic: Matrix) -> Expr | None
```

Internal implementation. Call `SemiMajorAxis(conic)` directly.

<a id="central_conic.SemiMinorAxis"></a>

## SemiMinorAxis Objects

```python
class SemiMinorAxis(Function)
```

Computes the semi-minor axis length of a conic.

The returned value is:
 - a positive number for ellipses;
 - infinity for parabolas;
 - an imaginary number for hyperbolas and complex ellipses;
 - nan for ideal point conics;
 - zero for the other degenerate conics.

Returns an unevaluated `sympy.Function` if we can't tell which axis is
shorter.

<a id="central_conic.SemiMinorAxis.eval"></a>

#### eval

```python
@classmethod
def eval(cls, conic: Matrix) -> Expr | None
```

Internal implementation. Call `SemiMinorAxis(conic)` directly.

<a id="circle"></a>

# circle

<a id="circle.Circle"></a>

#### Circle

```python
def Circle(center: Matrix | Sequence[Expr], radius: Expr) -> Matrix
```

Creates a circle from its center and radius.

<a id="circle.CircleRadius"></a>

#### CircleRadius

```python
def CircleRadius(circle: Matrix) -> Expr
```

Computes the radius of a circle conic.

The result is not specified if the conic matrix is not a circle.
The computation is based on
[research/director_circle.py](../src/research/director_circle.py).

<a id="circle.DirectorCircle"></a>

#### DirectorCircle

```python
def DirectorCircle(conic: Matrix) -> Matrix
```

Computes the director circle of a conic. It's also called orthoptic
circle or Fermat–Apollonius circle.

*Definition*: https://en.wikipedia.org/wiki/Director_circle<br>
*Formula*: [research/director_circle.py](../src/research/director_circle.py)

<a id="circle.UNIT_CIRCLE"></a>

#### UNIT\_CIRCLE

The circle at the origin with radius 1.

<a id="circle.COMPLEX_UNIT_CIRCLE"></a>

#### COMPLEX\_UNIT\_CIRCLE

The circle at the origin with radius 𝑖.

<a id="line"></a>

# line

<a id="line.LineContainsPoint"></a>

#### LineContainsPoint

```python
def LineContainsPoint(line: Matrix,
                      point: Matrix | Sequence[Expr]) -> bool | None
```

Tells whether `point` is on `line`.

Returns None if undecidable.

<a id="line.HorizontalLine"></a>

#### HorizontalLine

```python
def HorizontalLine(y: Expr) -> Matrix
```

Constructs a horizontal line with the given y-coordinate.

<a id="line.VerticalLine"></a>

#### VerticalLine

```python
def VerticalLine(x: Expr) -> Matrix
```

Constructs a vertical line with the given x-coordinate.

<a id="line.LineBetween"></a>

#### LineBetween

```python
def LineBetween(point1: Matrix | Sequence[Expr],
                point2: Matrix | Sequence[Expr]) -> Matrix
```

Connects two projective points with a line.

<a id="line.ParallelLine"></a>

#### ParallelLine

```python
def ParallelLine(with_line: Matrix,
                 through_point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a line through a point parallel to a line.

<a id="line.PerpendicularLine"></a>

#### PerpendicularLine

```python
def PerpendicularLine(to_line: Matrix,
                      through_point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a line through a point perpendicular to a line.

<a id="line.LineThroughPoint"></a>

#### LineThroughPoint

```python
def LineThroughPoint(point: Matrix | Sequence[Expr],
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

<a id="line.AngleBisector"></a>

#### AngleBisector

```python
def AngleBisector(line1: Matrix, line2: Matrix) -> Matrix
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

<a id="line.PerpendicularBisector"></a>

#### PerpendicularBisector

```python
def PerpendicularBisector(point1: Matrix | Sequence[Expr],
                          point2: Matrix | Sequence[Expr]) -> Matrix
```

Constructs the perpendicular bisector of two points.

<a id="line.AreParallel"></a>

#### AreParallel

```python
def AreParallel(line1: Matrix, line2: Matrix) -> bool | None
```

Tells whether line1 and line2 are parallel.

Returns True if they are parallel, False if not, None if undecidable.
Considers the ideal line parallel to everything.

<a id="line.ArePerpendicular"></a>

#### ArePerpendicular

```python
def ArePerpendicular(line1: Matrix, line2: Matrix) -> bool | None
```

Tells whether line1 and line2 are perpendicular.

Returns True if they are perpendicular, False if not, None if undecidable.
Considers the ideal line perpendicular to everything.

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

<a id="conic.ConicFromPoly"></a>

#### ConicFromPoly

```python
def ConicFromPoly(poly: Expr | Poly,
                  *,
                  x: Symbol = abc.x,
                  y: Symbol = abc.y) -> Matrix
```

Constructs the 3×3 symmetric matrix representation of a conic section
from a two-variable quadratic polynomial in the form of
```
ax² + bxy + cy² + dx + ey + f
```

The resulting matrix is:
```
[a, b/2, d/2]
[b/2, c, e/2]
[d/2, e/2, f]
```

<a id="conic.ConicThroughPoints"></a>

#### ConicThroughPoints

```python
def ConicThroughPoints(p1: Matrix | Sequence[Expr],
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

<a id="conic.ConicFromFocusAndDirectrix"></a>

#### ConicFromFocusAndDirectrix

```python
def ConicFromFocusAndDirectrix(focus: Matrix | Sequence[Expr],
                               directrix: Matrix,
                               eccentricity: Expr) -> Matrix
```

Constructs a conic from its focus, directrix and eccentricity.

*Formula*:
[research/conic_from_focus_and_directrix.py](../src/research/conic_from_focus_and_directrix.py)

<a id="conic.Eccentricity"></a>

#### Eccentricity

```python
def Eccentricity(conic: Matrix) -> Expr
```

Computes the eccentricity of a conic section.

The result is ambiguous in case of degenerate conics: evaluate
`(Eccentricity(conic), Eccentricity(-conic))` to get both values.

*Formula*: https://en.wikipedia.org/wiki/Conic_section#Eccentricity_in_terms_of_coefficients<br>
*Own research*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

<a id="conic.AxisDirection"></a>

#### AxisDirection

```python
def AxisDirection(conic: Matrix) -> Matrix
```

Returns the direction of the major axis of the conic section,
namely the ideal point in that direction.

Returns `[0, 0, 0]ᵀ` in case of circles.

The result is ambiguous in case of degenerate conics: evaluate
`(AxisDirection(conic), AxisDirection(-conic))` to get both values.

*Formula*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

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

<a id="conic.ProjectiveConicCenter"></a>

#### ProjectiveConicCenter

```python
def ProjectiveConicCenter(conic: Matrix) -> Matrix
```

Computes the generalized projective center of a conic.

It's equivalent to [ConicCenter](#central_conic.ConicCenter) (returns a
finite point) for
 - real and complex ellipses
 - hyperbolas
 - conics consisting of a single finite (Euclidean) point
 - crossing finite line pairs

For parabolas returns the ideal point on it
([proof](../src/research/parabola_center.py))

For other line pair conics and ideal point conics returns `[0, 0, 0]ᵀ`.

<a id="conic.PolePoint"></a>

#### PolePoint

```python
def PolePoint(conic: Matrix, polar_line: Matrix) -> Matrix
```

Computes the pole point of a conic with respect to the given polar line.

If the conic is degenerate, i.e. it factors into `l₁` and `l₂` real or
complex conjugate lines, the pole is
 - `[0, 0, 0]ᵀ` if `l₁`, `l₂`, and `polar_line` are concurrent;
 - the intersection of `l₁` and `l₂` otherwise.

*Pole / polar identity*: `conic * pole_point = polar_line`<br>
*Source*: https://en.wikipedia.org/wiki/Pole_and_polar#Calculating_the_pole_of_a_line

<a id="conic.PolarLine"></a>

#### PolarLine

```python
def PolarLine(conic: Matrix, pole_point: Matrix | Sequence[Expr]) -> Matrix
```

Computes the polar line of a conic with respect to the given pole point.

*Pole / polar identity*: `conic * pole_point = polar_line`
*Source*: https://en.wikipedia.org/wiki/Pole_and_polar#Calculating_the_pole_of_a_line

<a id="conic.ConicContainsPoint"></a>

#### ConicContainsPoint

```python
def ConicContainsPoint(conic: Matrix,
                       point: Matrix | Sequence[Expr]) -> bool | None
```

Checks if a point lies on a conic.

Returns None if undecidable.

<a id="conic.ConicContainsLine"></a>

#### ConicContainsLine

```python
def ConicContainsLine(conic: Matrix, line: Matrix) -> bool | None
```

Checks if a line lies on a conic.

Returns None if undecidable.

*Formula*:
[research/conic_line_containment.py](../src/research/conic_line_containment.py)

<a id="degenerate_conic"></a>

# degenerate\_conic

<a id="degenerate_conic.LinePair"></a>

#### LinePair

```python
def LinePair(line1: Matrix, line2: Matrix) -> Matrix
```

Constructs a conic section from two projective lines.

<a id="degenerate_conic.PointConic"></a>

#### PointConic

```python
def PointConic(point: Matrix | Sequence[Expr]) -> Matrix
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

Returns `[0, 0, 0]ᵀ` for double line pairs, or an unspecified 3d
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

<a id="distance"></a>

# distance

<a id="distance.PointPointDistance"></a>

#### PointPointDistance

```python
def PointPointDistance(point1: Matrix | Sequence[Expr],
                       point2: Matrix | Sequence[Expr]) -> Expr
```

Computes the signed distance between two points.

Special cases:
 - the distance between finite and ideal points is infinity
 - the distance between two ideal points is `nan`

<a id="distance.PointLineDistance"></a>

#### PointLineDistance

```python
def PointLineDistance(point: Matrix | Sequence[Expr], line: Matrix) -> Expr
```

Computes the signed distance between a point and a line.

Special cases:
 - the distance between finite points and the ideal line is infinity
 - the distance between ideal points and any projective lines is `nan`

<a id="sympy_utils"></a>

# sympy\_utils

<a id="sympy_utils.AddEq"></a>

#### AddEq

```python
def AddEq(*eqs: Eq) -> Eq
```

Adds multiple sympy equations.

<a id="sympy_utils.SubEq"></a>

#### SubEq

```python
def SubEq(eq0: Eq, eq1: Eq) -> Eq
```

Subtracts one sympy equation from another.

<a id="sympy_utils.MulEq"></a>

#### MulEq

```python
def MulEq(eq: Eq, factor: Expr | Eq) -> Eq
```

Multiplies a sympy equation by a factor or another equation.

<a id="sympy_utils.DivEq"></a>

#### DivEq

```python
def DivEq(eq: Eq, denom: Expr) -> Eq
```

Divides a sympy equation by a denominator.

<a id="sympy_utils.SwapEq"></a>

#### SwapEq

```python
def SwapEq(eq: Eq) -> Eq
```

Swaps the lhs and rhs of a sympy equation.

<a id="sympy_utils.EqChain"></a>

#### EqChain

```python
def EqChain(*expressions: Expr) -> Eq | Expr
```

Creates a chain of equations, i.e. `expr_1 = expr_2 = ... = expr_n`.

<a id="sympy_utils.FactorRadicals"></a>

#### FactorRadicals

```python
def FactorRadicals(expr: Expr) -> Expr
```

Factors all `Pow` and `sqrt` subexpressions inside `expr`.

<a id="conic_classification"></a>

# conic\_classification

<a id="conic_classification.IsDegenerate"></a>

#### IsDegenerate

```python
def IsDegenerate(conic: Matrix) -> bool | None
```

Tells whether the conic is degenerate.

Degenerate conics consist of a single projective point or a pair of
projective lines. The zero matrix is also considered degenerate.
Returns None if undecidable.

<a id="conic_classification.IsNonDegenerate"></a>

#### IsNonDegenerate

```python
def IsNonDegenerate(conic: Matrix) -> bool | None
```

Tells whether the conic is non-degenerate.

Non-degenerate conics include real or complex ellipses, parabolas and
hyperbolas. Returns None if undecidable.

<a id="conic_classification.IsFiniteConic"></a>

#### IsFiniteConic

```python
def IsFiniteConic(conic: Matrix) -> bool | None
```

Tells whether all points on the conic are finite.

Returns None if undecidable.

<a id="conic_classification.IsComplexEllipse"></a>

#### IsComplexEllipse

```python
def IsComplexEllipse(conic: Matrix) -> bool | None
```

Tells Whether the conic is an ellipse with a real center and imaginary
radii.

Returns None if undecidable.

<a id="conic_classification.IsEllipse"></a>

#### IsEllipse

```python
def IsEllipse(conic: Matrix) -> bool | None
```

Tells whether the conic is an ellipse with real radii.

Returns None if undecidable.

<a id="conic_classification.IsCircle"></a>

#### IsCircle

```python
def IsCircle(conic: Matrix) -> bool | None
```

Tells whether the conic is a circle.

Returns None if undecidable.

<a id="conic_classification.IsParabola"></a>

#### IsParabola

```python
def IsParabola(conic: Matrix) -> bool | None
```

Tells whether the conic is a parabola.

Returns None if undecidable.

<a id="conic_classification.IsHyperbola"></a>

#### IsHyperbola

```python
def IsHyperbola(conic: Matrix) -> bool | None
```

Tells whether the conic is a hyperbola.

Returns None if undecidable.

<a id="conic_classification.IsCircular"></a>

#### IsCircular

```python
def IsCircular(conic: Matrix) -> bool | None
```

Tells whether there is a single center point around which the conic is
invariant under all rotations.

Circles, complex circles, zero-radius circles have such circular symmetry.
Double ideal lines are not considered circular. Returns None if undecidable.

<a id="conic_classification.IsLinePair"></a>

#### IsLinePair

```python
def IsLinePair(conic: Matrix) -> bool | None
```

Tells whether the conic is the union of two projective lines.

Returns None if undecidable.

<a id="conic_classification.IsDoubleLine"></a>

#### IsDoubleLine

```python
def IsDoubleLine(conic: Matrix) -> bool | None
```

Tells whether the conic consists of two coincident projective lines.

Returns None if undecidable.

<a id="conic_classification.IsPointConic"></a>

#### IsPointConic

```python
def IsPointConic(conic: Matrix) -> bool | None
```

Tells whether the conic consists of a single projective point.

Returns None if undecidable.

A conic is a point conic iff it's degenerate and splits to two lines with
complex coordinates.

<a id="conic_classification.IsFinitePointConic"></a>

#### IsFinitePointConic

```python
def IsFinitePointConic(conic: Matrix) -> bool | None
```

Tells whether the conic consists of a single finite (Euclidean) point.

Returns None if undecidable.

<a id="point"></a>

# point

<a id="point.ORIGIN"></a>

#### ORIGIN

The point at (0, 0)

<a id="point.IdealPoint"></a>

#### IdealPoint

```python
def IdealPoint(x: Expr, y: Expr) -> Matrix
```

Creates an ideal point at the given direction.

<a id="point.IdealPointOnLine"></a>

#### IdealPointOnLine

```python
def IdealPointOnLine(line: Matrix) -> Matrix
```

Returns the coordinates of the ideal point on the line.

The first two coordinates specify the line's direction. The third one is
always zero.

If the line is the ideal line, returns `[0, 0, 0]ᵀ`.

<a id="point.PointToXY"></a>

#### PointToXY

```python
def PointToXY(point: Matrix | Sequence[Expr]) -> tuple[Expr, Expr]
```

Computes the Euclidean coordinates of a projective point.

<a id="point.PointToVec3"></a>

#### PointToVec3

```python
def PointToVec3(point: Matrix | Sequence[Expr]) -> Matrix
```

Computes the homogeneous coordinates of a projective point.

<a id="point.Centroid"></a>

#### Centroid

```python
def Centroid(*points: Sequence[Expr]) -> tuple[Expr, Expr]
```

Computes the centroid of a set of points.

<a id="point.PerpendicularFoot"></a>

#### PerpendicularFoot

```python
def PerpendicularFoot(point: Matrix | Sequence[Expr],
                      line: Matrix) -> tuple[Expr, Expr]
```

Computes the foot of the perpendicular through `point` to `line`.

Returns `(nan, nan)` when `point`, `line` or both are infinite.

*Formula*: https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

<a id="transform"></a>

# transform

<a id="transform.TransformConic"></a>

#### TransformConic

```python
def TransformConic(conic: Matrix, transformation: Matrix) -> Matrix
```

Applies a projective transformation on a conic.

<a id="transform.TransformLine"></a>

#### TransformLine

```python
def TransformLine(line: Matrix, transformation: Matrix) -> Matrix
```

Applies a projective transformation on a projective line.

<a id="transform.Translate"></a>

#### Translate

```python
def Translate(dx: Expr, dy: Expr) -> Matrix
```

Computes the transformation matrix for a 2d translation.

<a id="transform.Rotate"></a>

#### Rotate

```python
def Rotate(angle: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix
```

Computes the transformation matrix for a rotation around a point.

<a id="transform.ScaleXY"></a>

#### ScaleXY

```python
def ScaleXY(scale_x: Expr,
            scale_y: Expr,
            x0: Expr = 0,
            y0: Expr = 0) -> Matrix
```

Computes the projective transformation matrix for scaling along the x-
and y-axes.

<a id="transform.Scale"></a>

#### Scale

```python
def Scale(scale: Expr, x0: Expr = 0, y0: Expr = 0) -> Matrix
```

Computes the projective transformation matrix for a uniform scaling
transformation.

<a id="parabola"></a>

# parabola

<a id="parabola.ParabolaDirectrix"></a>

#### ParabolaDirectrix

```python
def ParabolaDirectrix(parabola: Matrix) -> Matrix
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

<a id="parabola.ParabolaFocus"></a>

#### ParabolaFocus

```python
def ParabolaFocus(parabola: Matrix) -> Matrix
```

Computes the focus of a parabola represented as a conic matrix.

Special cases for other conic types:
- Returns `[0, 0, 0]ᵀ` for
  - conics with one ideal point;
  - conics containing the ideal line.
- Raises `ValueError` if the conic provably has 0 or 2 ideal points.
- Returns an unspecified 3D column vector in all other cases.

*Formula*:
[research/focus_directrix_eccentricity.py](../src/research/focus_directrix_eccentricity.py)

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

<a id="polar_conic.PointAtAngle"></a>

#### PointAtAngle

```python
def PointAtAngle(polar: Matrix, theta: Expr) -> Matrix
```

Computes the coordinates of the projective point on a polar conic
corresponding to a certain angle.

<a id="polar_conic.ConicFromPolarMatrix"></a>

#### ConicFromPolarMatrix

```python
def ConicFromPolarMatrix(polar: Matrix) -> Matrix
```

Transforms a conic from polar to quadratic form.

The algorithm is essentially applying the polar matrix as a projective
transformation on the unit circle.

<a id="ellipse"></a>

# ellipse

<a id="ellipse.Ellipse"></a>

#### Ellipse

```python
def Ellipse(center: Matrix | Sequence[Expr],
            r1: Expr,
            r2: Expr,
            *,
            r1_angle: Expr = None,
            r1_direction: Expr = None) -> Matrix
```

Constructs an ellipse from its center, radii and the rotation angle of
the first radius.

*Formula*:
[research/ellipse_from_params.py](../src/research/ellipse_from_params.py)

<a id="ellipse.EllipseFromFociAndPoint"></a>

#### EllipseFromFociAndPoint

```python
def EllipseFromFociAndPoint(focus1: Matrix | Sequence[Expr],
                            focus2: Matrix | Sequence[Expr],
                            point: Matrix | Sequence[Expr]) -> Matrix
```

Constructs an ellipse from its focus points and an incident point.

<a id="ellipse.SteinerEllipse"></a>

#### SteinerEllipse

```python
def SteinerEllipse(point1: Matrix | Sequence[Expr],
                   point2: Matrix | Sequence[Expr],
                   point3: Matrix | Sequence[Expr]) -> Matrix
```

Constructs the Steiner circumellipse for the given points.

The ellipse goes through the three points and is centered at the triangle's
centroid.

*Definition*: https://en.wikipedia.org/wiki/Steiner_ellipse<br>
*Formula*: [research/steiner_ellipse.py](../src/research/steiner_ellipse.py)

<a id="ellipse.SteinerInellipse"></a>

#### SteinerInellipse

```python
def SteinerInellipse(point1: Matrix | Sequence[Expr],
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

<a id="intersection.LineXLine"></a>

#### LineXLine

```python
def LineXLine(line1: Matrix, line2: Matrix) -> Matrix
```

Computes the intersection of two lines.

Returns an ideal point if the lines are parallel, or `[0, 0, 0]ᵀ` if they
coincide.

<a id="intersection.ConicXLine"></a>

#### ConicXLine

```python
def ConicXLine(
    conic: Matrix, line: Matrix
) -> tuple[Matrix | Sequence[Expr], Matrix | Sequence[Expr]] | NaN
```

Intersects a conic with a line. Returns two points.

Special cases:
 - The intersection points coincide if the line is tangent to the conic.
 - They are complex conjugates if the line doesn't intersect the conic at
   a real point.
 - Returns an unevaluated `sympy.Function` for symbolic conics.
 - Returns None if the conic contains the entire line.

*Algorithm*: Jürgen Richter-Gebert, Perspectives on Projective Geometry,
section 11.3

