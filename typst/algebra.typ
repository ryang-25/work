#import "notes.typ": *
#show: great-theorems-init

#show: notes.with( title: "Algebra", subject: "Notes from the book", author: "Roland Yang")
#show heading: reset-counter(tcounter, levels: 2)

#let end = math.op("End")

#heading(numbering: none, "Introduction")

We should start with Chapters 1-5, then Chapter 8, then Chapter 11. Try all the
exercises in the book.

= The Integers

#problem[Find an example illustrating why the hypothesis that $a != 0$ is
 necessary in the statement of Lemma 1.2.]
#solution[We have a trivial example. If $a = 0$ and $b = 1$, then $b divides a$,
 as $a = 0b$.]

#problem[Let $n$ be a positive integer, and let $a$ and $b$ be integers. Prove
 that $a$ and $b$ have the same remainder when divided by $n$ if and only if
 $a - b = n k$ for some integer $k$.]

= Modular Arithmetic

= Rings
= The Category of Rings
= Canonical Decomposition, Quotients, and Isomorphism Theorems

= Integral Domains

#problem[]
#problem[]
#problem[]
#problem[]
#problem[]
#problem[Prove that if $R$ is an integral domain, then $(x)$ is a prime ideal in
 $R[x]$.]
#solution[By the first isomorphism theorem, $R[x]\/(x) tilde.equiv R$, $R[x]\/
 (x)$ is an integral domain, and $(x)$ is a prime ideal of $R[x]$.]
#problem[]
#problem[]
#problem[Prove that the ideal $(y+1,x^2+1)$ is maximal in $RR[x,y]$ and contains
 the ideal $(y-x^2)$. This ideal does not correspond to a point with real
 coordinates on the parabola $y-x^2$. (This phenomenon is one reason why
 algebraic geometry is ‘easier’ over the complex numbers.)]

#counter(heading).update(7)

= Modules and Abelian Groups

- A vector space over a field $k$ is a set $V$ endowed with $(+, dot)$ by
  scalars
- The polynomial ring $k[x]$ over $k$ is a $k$-vector space
- This includes the ring addition axioms as well as:
  - $forall v in V, 1_k v = v$
  - $forall a, b in k, forall v in V, (a b)v = a (b V)$
  - $ forall a, b in k, forall v, w in V, (a + b)v  &= a v + b v \ a( v + w) &=
    a v + a w $
- Scalar multiplication defines an _action_ of $k$ on $V$ (since these are not
  really the ring axioms, but they do look very similar)
- We see ideal absorption $(forall r in R)(forall a in I), a r in I$ is similar
  to scalar multiplication:\ $R times I arrow.r I$
- Think of ideal over ring as some vector space
  - We define an $R$-module as this generalization, notated $M$
- A vector space is a "$k$-module" (my words)—an $R$-module where $R$ is a field
  $k$
- Rings are (obviously) $R$-modules
- Define a ring homomorphism $f: R arrow.r S$. Define an action of $R$ on $S$ by
  $forall r in R, forall s in S, r dot s = f(r)s$.
  - This is an $R$-module over $S$
- $R\/I$ is also an $R$-module. $r dot (1 + I)$ generates $R\/I$
  - Then this is a cyclic $R$-module
- If ground ring is $ZZ$ then we can have the action $n dot a = underbrace(
  (1+ dots.c + 1), n "times") dot a = underbrace(1 dot a + dots.c + 1 dot a,
  n "times")$

#heading(numbering: none, level: 2, "Exercises")
#problem[ Recall that the set of _linear transformations_ $V arrow V$ of a
 $k$-vector space $V$ forms a ring, call this ring $end(V)$. Prove that
 multiplication by scalars determines a ring homomorphism $k arrow end(V)$. ]
#solution[
  It is simple to verify the definition of a homomorphism. For every $c in k$,
  we define \ $forall v in V$, $mu(v) = c v$. Then we see

  $
    mu(c_1 + c_2) &= (c_1 + c_2)v = c_1v + c_2 v = mu(c_1) + mu(c_2) \
    mu(c_1 c_2) &= (c_1c_2)v = c_1(c_2 v) = mu(c_1mu(c_2)) = mu(c_1) compose mu(c_2)
  $
  and the identity $1_k arrow.r.long.bar id_V$.
  
]

#problem[ Verify that if $f : R arrow S$ is a ring homomorphism, then $S$ is an
 $R$-module in a natural way. ]
#solution[
  From the action defined in Example 8.6, we can verify axioms (i)-(v) by virtue of the ring homomorphism, then it suffices to show that for an $a, b in S$
  $
    (r s) dot a = f(r s) dot a = f(r)f(s)a = f(r)(f(s) a) = r dot (s dot a) \
    (r + s) dot a = f(r + s) = (f(r) + f(s))a  = f(r)a + f(s)a = r dot a + s dot a \
    r dot (a + b) = f(r)(a + b) = f(r)a + f(r)b =   r dot a + r dot b
  $
  by the ring axioms in Definition 3.1.
]

#problem[ Let $K$ be a field, and let $k subset.eq K$ be a subfield. Show that
 $K$ is a $k$-vector space in a natural way. ]
#solution[ We have the inclusion $iota : k arrow K$, so we can view $K$ as a
 $k$-module by exercise 2, and this is simply a $k$-vector space. ]

#problem[
  Convince yourself of the truth of Proposition 8.10.
]


#problem[ We have seen that every abelian group, a.k.a. $ZZ$-module, has
 a _unique_ $ZZ$-module structure. For example, $ZZ$ itself can be realized as
 a $ZZ$-module in only one way. Prove that $ZZ$ can be realized as a $ZZ
 [x]$ module in _infinitely many_ different ways. ]
#solution[]

#problem[ Prove that if $R$ is a ring, $r in R$ is in the center of $R$ (i.e., it
 commutes with every element of $R$), and $M$ is an $R$-module, then the
 function $M arrow M$ defined by $m arrow.bar r m$ is an $R$-module
 homomorphism. Find an example showing that something may go wrong if $r$ is
 not in the center of $R$. ]
#solution[ Conversely, if $r$ is not in the center, then we cannot say that $n r
 a = r n a$, and indeed we can see this with most matrix multiplications. ]

#problem[
  Prove that the composition of two $R$-module homomorphisms is an $R$-module
  homomorphism.
]
#solution[
  Let $f : M arrow N$ and $g : N arrow O$ be $R$-module homomorphisms. Then
  composing is verifying that
  $
    (g compose f)(a+b) &= g(f(a + b)) = g(f(a) + f(b)) = g(f(a)) + g(f(b)) =
    (g compose f )(a) + (g compose f)(b) \
    (g compose f)(a dot b) &= g(f(a dot b)) = g(f(a) dot f(b)) =
    g(f(a)) dot g(f(b)) = (g compose f )(a) dot (g compose f)(b) \
  $
]

#problem[
  Let $f : M arrow N$ be a homomorphism of $R$-modules. Prove that $f(0_M) = f(0_N).$
]
#solution[
  This is exactly what we did for rings.
  $
    f(0_M) = f(0_M + 0_M) = f(0_M) + f(0_M)
  $
  and then $f(0_M) = 0_N$ by cancellation.
]

#problem[
  Does the category of $R$-modules have initial objects? Final objects?
]
#solution[]

#tcounter.update(13)
#problem[
  Let $R$ be a ring. A nonzero $R$-module $M$ is 'simple' if its only submodules
  are ${0}$ and $M$. Let $M,N$ be simple modules, and let $phi.alt : M arrow N$
  be a homomorphism of $R$-modules. Prove that either $phi.alt = 0$ or $phi.alt$
  is an isomorphism. (This statement is known as _Schur's lemma._)
]

#counter(heading).update(10)
= Groups—Preliminaries
