#import "notes.typ": *
#show: great-theorems-init

#show: notes.with(
  title: "Ideals, Varieties, and Algorithms",
  subject: "Exercises and Solutions",
  author: "Roland Yang",
)
#show heading: reset-counter(tcounter, levels: 2)

#let upbold(x) = $bold(upright(#x))$
#let I = $I$
#let k
#let kx = $k[x]$
#let ky = $k[y]$
#let kxy = $k[x,y]$
#let kpolyx = $k[x_1, dots, x_n]$
#let kpolyy = $k[y_1, dots, y_m]$
#let fs = $f_1, dots, f_m$

/*
#solution[
  #lifblock[For a $J subset.eq k[V]$, $V_v (J) = emptyset$, let
   $accent(J, tilde)$ be the corresponding ideal in $k[x_1,dots,x_n]$. By
   definition, $V_v (J) = bold(V)(accent(J, tilde)) = emptyset$, which by the
   Weak Nullstellensatz implies $accent(J, tilde) = k[x_1,dots,x_n]$. By the
   isomorphism, $J = k[V]$.]

  #rifblock[$J = k[V]$ implies
    $accent(J, tilde) = k[x_1,dots,x_n]$ and $bold(V) (I) = emptyset$ for an
    infinite ring. (I couldn't think of another proof.
  ]
]
*/

// Exercise 2 (using the fundamental theorem) we have that when $J = angle.l phi.alt angle.r $ that $phi.alt$ is a nonzero constant, so the conclusion is the same.

= Geometry, Algebra, and Algorithms
= Gröbner Bases
= Elimination Theory
= The Algebra-Geometry Dictionary
= Polynomial and Rational Functions on a Variety

==
==
==
==
==
== Relative Finiteness and Noether Normalization

#problem[]

#tcounter.update(5)
#problem[
  Show that $kxy\/angle.l x^a - y^b angle.r$ is finite both over #kx and #ky,
  but that $kxy\/angle.l x^(a+1) - x y^b angle.r$ is finite over #ky but not
  #kx. Interpret this geometrically when #k is algebraically closed.
]
#solution[
  We can see that both
  $kx inter kxy\/angle.l x^a - y^b angle.r = ky inter kxy\/angle.l x^a - y^b angle.r = {0}$
  and every element can be expressed both as a #kx\-linear combo in
  $1,dots,y^(b-1)$ and a #ky\-linear combo in $1,dots,x^(a-1)$.
]

#problem[

]


#counter(heading).step()
= Invariant Theory of Finite Groups

==

#problem[]

#solution[

  #h(0.5em) $arrow.long$ By Definition 1, if $f$ symmetric then $f(x,y,z) = f(y,x,z) = f(y,z,x)$. \
  #h(0.5em) $arrow.l.long$ If $y$ and $x$ are interchangeable and $x$ and $z$ are
  interchangeable then $y$ and $z$ are transitively interchangeable, and with these
  you can generate all $3!$ permutations.
]

#tcounter.step()

#problem[]
#problem[]

#tcounter.update(15)

#problem[]
#solution[Since $x y$ is of degree 2 it must be a product of degree 1 polynomials.
  Then by definition the power sum $s_1 = x + y$, and
  $s_1^2 = x^2 + 2 x y + y^2 = x^2 + y^2$, where the powers are larger and not
  divisible by $x y$, so it is impossible!]

// 7.2
== Finite Matrix Groups and Rings of Invariants

#problem[]

#problem[Suppose that $A in "GL"(n,k)$ satisfies $A^m = I_n$ for some positive integer $m$. ]
#solution[By contradiction, if there exists $A^a = A^b, a < b$, then $m$
  would not be the smallest integer $A^m = I_n$. Every power $p > m$ can
  be reduced modulo $m$, proving closure. ]


#problem[Write down the six permutation matrices of GL$(3, k)$.]
#solution[$
  mat(1, 0, 0; 0, 1, 0; 0, 0, 1),
  mat(1, 0, 0; 0, 0, 1; 0, 1, 0),
  mat(0, 1, 0; 1, 0, 0; 0, 0, 1),
  mat(0, 1, 0; 0, 0, 1; 1, 0, 0),
$]

#tcounter.step()
#problem[]

#tcounter.step()
#problem[]

#problem[
]
#solution[
  The dual of the tetrahedron is the tetrahedron (what other regular polyhedron has 4 vertices and 4 faces?)
]

#tcounter.update(13)
#problem[]
#solution[
  + $ mat(i, 0; 0, -i), mat(-1, 0; 0, -1), mat(-i, 0; 0, i), mat(1, 0; 0, 1) $
  + $f(x,y) = f(-i x, i y) = f(-x, -y) = f(i x, -y)$ implies that if
    $f = sum_(i j) a_(i j) x^i y^j$ then we see that for the imaginary components:
    $ (-1)^i i^(i+j) = (-1)^j i^(i+j) = 1 $ which implies that $i + j$ must either
    be $0$ or $2 mod 4$; for the first case, $i$ and $j$ must be even and
    necessarily $0$ or $2$ mod 4. For the latter case, $i$ or $j$ must be odd,
    implying that $i = j =$ $1$ or $3 mod 4$, so $i$ and $j$ must be congruent $mod 4$.
]

== Generators for the Ring of Invariants

#problem[]

// 7.4
== Relations Among Generators and the Geometry of Orbits

#problem[
  Given $fs in kpolyx$, let $I = {g in kpolyy | g(fs) = 0}$.
  + Prove that #I is an ideal of #kpolyy.
  + If $f in kpolyx$ and $f = g(fs)$ is one representation of $f$ in terms of
    #fs, provat that all such representations are given by $f = g(fs) + h(fs)$
    as $h$ varies over #I.
]
#solution[
  + To prove that #I is an ideal, we see that it is closed under addition
    ($h_1 + h_2 = 0 + 0 = 0$), the zero element is trivially included, and
    each $h$ is its own additive inverse where $h - h = 0 - 0 = 0$, and similarly
    absoprtion where $0 dot r = r dot 0 = 0$.
  + By contradiction, suppose $exists g'$ where $g' - g = r$, where
    $r(fs) eq.not 0$ in #kpolyx, but this is a contradiction since this implies
    $f - f eq.not 0$ in #kpolyx!
]

#problem[
  Let $fs in kpolyx$ and let $I subset.eq kpolyy$ be the ideal of relations defined
  in Exercise 1.
  + Prove that the map sending a coset $[g]$ to $g(fs)$ defines a well-defined ring
    homomorphism $ phi.alt : kpolyy\/I --> k[fs]. $
  + Prove that the map $phi.alt$ of part (a) is one-to-one and onto. Thus
    $phi.alt$ is a ring isomorphism.
]
#solution[
  +
  + It is one-to-one: $ker phi.alt = 0$ and onto. Every function in $k[fs]$ has a
    corresponding coset.
  + Since $Phi$ is surjective and $ker Phi$ is exactly $I$ as it is defined in
    (1), $kpolyy\/I tilde.equiv k[fs]$ by the first ring isomorphism theorem.
]

#tcounter.update(4)
#problem[
  Complete Example 5 by showing that $I_F subset.eq k[u,v,w]$ is given by
  $I_F = angle.l u^2 w - v^2 - 4w^2 angle.r$ when
  $F = (x^2 + y^2, x^3 y - x y^3, x^2 y^2)$.
]
#solution[
  Using Sage this fact can be trivially verified.
]

#problem[]

#tcounter.update(8)
#problem[

]
#solution[
  - _Reflexive:_ $upbold(a) = I_n upbold(a)$.
  - _Symmetric:_ $upbold(b) = A upbold(a) -> A^(-1) upbold(b) = I_n upbold(a)$.
  - _Transitive:_ If $upbold(b) = A upbold(a)$ and $upbold(c) = B upbold(b)$,
    then $upbold(c) = (B A)upbold(a)$.
]

#problem[]
#solution[
  + The 1-orbit is at the origin, the 6-orbit is through the center of a face,
    the 8-orbit is through a corner, the 12-orbit is through an edge, and the
    24 orbit is through a non-centered point on a face.
]

= Projective Algebraic Geometry

== The Projective Plane

#let PRP = $PP^2(RR)$
#let Hinf = $H_infinity$

#problem[
  Using #PRP as given in Definition 1, we saw that the projective lines in
  $PP^2(RR)$ are $overline(L) = L union [L]_infinity$, and the line at $infinity$.
  + Prove that _any_ two distinct points in $PRP$ determine a unique projective
    line. Hint: There are three cases, depending on how many of the points are points
    at $infinity$.
  + Prove that _any_ two distinct projective lines in #PRP meet at a
    unique point. Hint: Do this case-by-case.
]

#tcounter.update(3)
#problem[
  This exercise will study what the hyperbola $x^2 - y^2 = 1$ looks like at $infinity$.
  + Explain why the equation $x^2 - y^2 = z^2$ gives a well-defined curve $C$ in #PRP.
    Hint: See the discussion following Definition 3.
]

#problem[]

#problem[
  When we use the $(x,y)$ coordinate system inside #PRP, we only view a piece of
  the projective plane. In particular, we miss the line at $infinity$. As in the
  text, we can use $(x,z)$ coordinates to view the line at $infinity$. Show that
  there is exactly one point in #PRP that is visible in neither $(x,y)$ nor $(x,z)$
  coordinates. How can we view what is happening at this point?
]
#solution[
  For a point we have the one-to-one maps $(x,y) mapsto (x,y,1)$ and
  $(x,z) mapsto (x,1,z)$. A point is then invisible if it is of the form
  $(x,0,0)$ and using homogeneous coordinates, it is the single point $(1,0,0)$
  in #PRP.
]

#problem[
  In the proof of Proposition 4, show that the image of the map (1) is disjoint
  from #Hinf.
]
#solution[
  Since the map sets $z = 1$, and there is no point $(0:0:0) in PRP$, the two
  are disjoint.
]

#problem[
  As in the text, the line #Hinf is defined by $z = 0$. Thus, points on #Hinf have
  homogeneous coordinates $(a:b:0)$, where $(a,b) eq.not (0,0)$.
  + A vertical affine line $x = c$ gives the projective line $x = c z$. Show that
    this meets #Hinf at the point $(0:1:0)$.
  + Show that a point on #Hinf different from $(0:1:0)$ can be written uniquely in
    the form $(1:m:0)$ for some real number $m$.
]
#solution[
  + At $x = c z |_(z = 0), 0 = c(0)$. Since $y$ is free, under homogeneous
    coordinates it meets at $(0:1:0)$.
  + For any point $(a:b:0)$ where $a eq.not 0$, we trivially see it can be
    written as $(1:b/a:0)$.
]

== Projective Space and Projective Varieties
