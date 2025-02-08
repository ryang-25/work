#import "notes.typ": *
#show: great-theorems-init

#show: notes.with(
  title: "Ideals, Varieties, and Algorithms",
  subject: "Notes from the book",
  author: "Roland Yang",
)

#solution[

  #h(0.5em) $arrow.long$ For a $J subset.eq k[V]$, $V_v (J) = emptyset$, let
   $accent(J, tilde)$ be the corresponding ideal in $k[x_1,dots,x_n]$. By
   definition, $V_v (J) = bold(V)(accent(J, tilde)) = emptyset$, which by the
   weak nullstellensatz implies $accent(J, tilde) = k[x_1,dots,x_n]$. By the
   isomorphism, $J = k[V]$.
  #h(0.5em) $arrow.l.long$ $J = k[V]$ implies $accent(J, tilde) = k
   [x_1,dots,x_n]$ and $bold(V) (I) = emptyset$ for an infinite ring. (I couldn't
   think of a proof by contrapositive.) 
]

// Exercise 2 (using the fundamental theorem) we have that when $J = angle.l phi.alt angle.r $ that $phi.alt$ is a nonzero constant, so the conclusion is the same.

=
= 
=
=
=

#counter(heading).step()

= Invariant Theory of Finite Groups

==

#problem[]

#solution[

  #h(0.5em) $arrow.long$ By Definition 1, if $f$ symmetric then $f(x,y,z) = f(y,x,z) = f(y,z,x)$. \
  #h(0.5em) $arrow.l.long$ If $y$ and $x$ are interchangeable and $x$ and $z$ are interchangeable then $y$ and $z$ are transitively interchangeable, and with these you can generate all $3!$ permutations
]

#tcounter.step()

#problem[]
#problem[]

#tcounter.update(15)

#problem[]
#solution[Since $x y$ is of degree 2 it must be a product of degree 1 polynomials. Then by definition the power sum $s_1
 = x + y$, and $s_1^2 = x^2 + 2 x y + y^2 = x^2 + y^2$, where the powers are larger and not divisible by $x y$, so it
 is impossible!]


== Finite Matrix Groups and Rings of Invariants

#problem[]

#problem[Suppose that $A in "GL"(n,k)$ satisfies $A^m = I_n$ for some positive integer $m$. ]
#solution[ By the contrapositive, if there exists $A^a = A^b, a < b$, then $m$ would not be the smallest integer $A^m =
 I_n$. Every power $p > m$ can be reduced modulo $m$, proving closure. ]


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
  #set enum(numbering: "a)", indent: 0.5em)
  + $ mat(i, 0; 0, -i), mat(-1, 0; 0, -1), mat(-i, 0; 0, i), mat(1, 0; 0, 1) $
  + $f(x,y) = f(-i x, i y) = f(-x, -y) = f(i x, -y)$ implies that if $f = sum_(i j) a_(i j) x^i y^j$ then we see that
    for the imaginary components: $ (-1)^i i^(i+j) = (-1)^j i^(i+j) = 1 $ which implies that $i + j$ must either be $0$
    or $2 mod 4$; for the first case, $i$ and $j$ must be even and necessarily $0$ or $2$ mod 4. For the latter case,
    $i$ or $j$ must be odd, implying that $i = j =$ $1$ or $3 mod 4$, so $i$ and $j$ must be congruent $mod 4$. ]