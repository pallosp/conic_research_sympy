# Table of Contents

* [matrix](#matrix)
  * [IsNonZeroMultiple](#matrix.IsNonZeroMultiple)
  * [MaxEigenvalue](#matrix.MaxEigenvalue)
  * [MinEigenvalue](#matrix.MinEigenvalue)
  * [ConicMatrix](#matrix.ConicMatrix)
  * [QuadraticForm](#matrix.QuadraticForm)
  * [SkewMatrix](#matrix.SkewMatrix)
  * [NonZeroCross](#matrix.NonZeroCross)
  * [IsDefinite](#matrix.IsDefinite)
* [central\_conic](#central_conic)
  * [ConicFromCenterAndPoints](#central_conic.ConicFromCenterAndPoints)
  * [ConicCenter](#central_conic.ConicCenter)
  * [SemiAxisLengths](#central_conic.SemiAxisLengths)
  * [SemiMajorAxis](#central_conic.SemiMajorAxis)
  * [SemiMinorAxis](#central_conic.SemiMinorAxis)
* [circle](#circle)
  * [Circle](#circle.Circle)
  * [CircleRadius](#circle.CircleRadius)
  * [DirectorCircle](#circle.DirectorCircle)
* [line](#line)
  * [HorizontalLine](#line.HorizontalLine)
  * [VerticalLine](#line.VerticalLine)
  * [LineBetween](#line.LineBetween)
  * [ParallelLine](#line.ParallelLine)
  * [PerpendicularLine](#line.PerpendicularLine)
  * [PerpendicularBisector](#line.PerpendicularBisector)
  * [AreParallel](#line.AreParallel)
  * [ArePerpendicular](#line.ArePerpendicular)
* [conic](#conic)
  * [ConicFromPoly](#conic.ConicFromPoly)
  * [ConicThroughPoints](#conic.ConicThroughPoints)
  * [ConicFromFocusAndDirectrix](#conic.ConicFromFocusAndDirectrix)
  * [Eccentricity](#conic.Eccentricity)
  * [AxisDirection](#conic.AxisDirection)
  * [SplitToLines](#conic.SplitToLines)
  * [IdealPoints](#conic.IdealPoints)
* [degenerate\_conic](#degenerate_conic)
  * [LinePair](#degenerate_conic.LinePair)
* [distance](#distance)
  * [PointPointDistance](#distance.PointPointDistance)
  * [PointLineDistance](#distance.PointLineDistance)
* [conic\_classification](#conic_classification)
  * [IsDegenerate](#conic_classification.IsDegenerate)
  * [IsNonDegenerate](#conic_classification.IsNonDegenerate)
  * [IsFiniteConic](#conic_classification.IsFiniteConic)
  * [IsComplexEllipse](#conic_classification.IsComplexEllipse)
  * [IsEllipse](#conic_classification.IsEllipse)
  * [IsParabola](#conic_classification.IsParabola)
  * [IsHyperbola](#conic_classification.IsHyperbola)
  * [IsCircular](#conic_classification.IsCircular)
  * [IsLinePair](#conic_classification.IsLinePair)
  * [IsDoubleLine](#conic_classification.IsDoubleLine)
* [point](#point)
  * [IdealPoint](#point.IdealPoint)
  * [IdealPointOnLine](#point.IdealPointOnLine)
  * [PointToXY](#point.PointToXY)
  * [PointToVec3](#point.PointToVec3)
  * [Centroid](#point.Centroid)
* [polar\_conic](#polar_conic)
  * [PointAtAngle](#polar_conic.PointAtAngle)
  * [ConicFromPolarMatrix](#polar_conic.ConicFromPolarMatrix)
* [ellipse](#ellipse)
  * [Ellipse](#ellipse.Ellipse)
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
                      m2: Matrix | Sequence[Expr]) -> bool
```

Tells whether two matrices are non-zero scalar multiples of each other.

Treats lists and tuples as column vectors.

<a id="matrix.MaxEigenvalue"></a>

#### MaxEigenvalue

```python
def MaxEigenvalue(symmetric_matrix2x2: Matrix) -> Expr
```

Returns the higher eigenvalue of a 2x2 symmetric matrix.

<a id="matrix.MinEigenvalue"></a>

#### MinEigenvalue

```python
def MinEigenvalue(symmetric_matrix2x2: Matrix)
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

Formula: `vᵀ·M·v`

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

<a id="central_conic.ConicCenter"></a>

#### ConicCenter

```python
def ConicCenter(conic: Matrix) -> Tuple[Expr, Expr]
```

Computes the center point of a conic.

Formula: ChatGPT

<a id="central_conic.SemiAxisLengths"></a>

#### SemiAxisLengths

```python
def SemiAxisLengths(conic: Matrix) -> Tuple[Expr, Expr]
```

Computes the semi-axis lengths of a conic.

Formula: ChatGPT

<a id="central_conic.SemiMajorAxis"></a>

#### SemiMajorAxis

```python
def SemiMajorAxis(conic: Matrix) -> Expr
```

Computes the semi-major axis length i.e. the center-vertex distance of
a conic.

It's infinity for parabolas, zero for degenerate conics, and imaginary for
complex ellipses.

<a id="central_conic.SemiMinorAxis"></a>

#### SemiMinorAxis

```python
def SemiMinorAxis(conic: Matrix) -> Expr
```

Computes the semi-minor axis length of a conic.

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
The computation is based on `research/director_circle.py`.

<a id="circle.DirectorCircle"></a>

#### DirectorCircle

```python
def DirectorCircle(conic: Matrix) -> Matrix
```

Computes the director circle of a conic. It's also called orthoptic
circle or Fermat–Apollonius circle.

Definition: https://en.wikipedia.org/wiki/Director_circle<br>
Formula: `research/director_circle.py`

<a id="line"></a>

# line

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
def ParallelLine(withLine: Matrix,
                 throughPoint: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a line through a point parallel to a line.

<a id="line.PerpendicularLine"></a>

#### PerpendicularLine

```python
def PerpendicularLine(toLine: Matrix,
                      throughPoint: Matrix | Sequence[Expr]) -> Matrix
```

Constructs a line through a point perpendicular to a line.

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

<a id="conic"></a>

# conic

<a id="conic.ConicFromPoly"></a>

#### ConicFromPoly

```python
def ConicFromPoly(poly: Expr | Poly,
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
def ConicThroughPoints(p1, p2, p3, p4, p5) -> Matrix
```

Computes the conic that goes through the given points.

Returns the conic matrix, or a zero matrix if the result is ambiguous.
The result is unique if no 2 points coincide and no 4 points are collinear.
The result is non-degenerate if no 3 points are collinear.

Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 10.1

<a id="conic.ConicFromFocusAndDirectrix"></a>

#### ConicFromFocusAndDirectrix

```python
def ConicFromFocusAndDirectrix(focus: Matrix, directrix: Matrix,
                               eccentricity: Expr) -> Matrix
```

Constructs a conic from its focus, directrix and eccentricity.

Formula: `research/conic_from_focus_and_directrix.py`

<a id="conic.Eccentricity"></a>

#### Eccentricity

```python
def Eccentricity(conic: Matrix)
```

Computes the eccentricity of a conic section.

The result is ambiguous in case of degenerate conics: evaluate
`(Eccentricity(conic), Eccentricity(-conic))` to get both values.

Formula: https://en.wikipedia.org/wiki/Conic_section#Eccentricity_in_terms_of_coefficients
Own research: `research/focus_directrix_eccentricity.py`

<a id="conic.AxisDirection"></a>

#### AxisDirection

```python
def AxisDirection(conic: Matrix) -> Matrix
```

Returns the direction of the major axis of the conic section,
namely the ideal point in that direction.

The result is ambiguous in case of degenerate conics: evaluate
`(AxisDirection(conic), AxisDirection(-conic))` to get both values.

Formula: `research/focus_directrix_eccentricity.py`

<a id="conic.SplitToLines"></a>

## SplitToLines Objects

```python
class SplitToLines(Function)
```

Splits a degenerate conic into two lines.

Special cases:
 - For non-degenerate conics the result is unspecified.
 - For point conics the lines will be complex conjugates.
 - For symbolic conics returns an unevaluated `sympy.Function`.

Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 11.1

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

<a id="degenerate_conic"></a>

# degenerate\_conic

<a id="degenerate_conic.LinePair"></a>

#### LinePair

```python
def LinePair(line1: Matrix, line2: Matrix) -> Matrix
```

Constructs a conic section from two projective lines.

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
 - the distance between Euclidean and ideal points is infinity
 - the distance between two ideal points is `nan`

<a id="distance.PointLineDistance"></a>

#### PointLineDistance

```python
def PointLineDistance(point: Matrix | Sequence[Expr], line: Matrix) -> Expr
```

Computes the signed distance between a point and a line.

Special cases:
 - the distance between Euclidean points and an ideal lines is infinity
 - the distance between ideal points and an any projective lines is `nan`

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

<a id="point"></a>

# point

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

If the line is the ideal line, returns a zero vector.

<a id="point.PointToXY"></a>

#### PointToXY

```python
def PointToXY(point: Matrix | Sequence[Expr]) -> Tuple[Expr, Expr]
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
def Centroid(*points: Sequence[Expr]) -> Tuple[Expr, Expr]
```

Computes the centroid of a set of points.

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

Formula: `research/ellipse_from_params.py`

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

Definition: https://en.wikipedia.org/wiki/Steiner_ellipse<br>
Formula: `research/steiner_ellipse.py`

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

Definition: https://en.wikipedia.org/wiki/Steiner_inellipse<br>
Formula: `research/steiner_ellipse.py`

<a id="intersection"></a>

# intersection

<a id="intersection.LineXLine"></a>

#### LineXLine

```python
def LineXLine(line1: Matrix, line2: Matrix) -> Matrix
```

Computes the intersection of two lines.

Returns an ideal point if the lines are parallel, or a zero vector if they
coincide.

<a id="intersection.ConicXLine"></a>

#### ConicXLine

```python
def ConicXLine(
    conic: Matrix, line: Matrix
) -> Tuple[Matrix | Sequence[Expr], Matrix | Sequence[Expr]] | NaN
```

Intersects a conic with a line. Returns two points.

Special cases:
 - The intersection points coincide if the line is tangent to the conic.
 - They are complex conjugates if the line doesn't intersect the conic at
   a real point.
 - Returns an unevaluated `sympy.Function` for symbolic conics.
 - Returns None if the conic contains the entire line.

Algorithm: Jürgen Richter-Gebert, Projective Geometry, section 11.3

