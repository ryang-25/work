#import "notes.typ": *
#show: great-theorems-init

#show: notes.with(
  title: "Algebra",
  subject: "Notes from the book",
  author: "Roland Yang",
)

#heading(numbering: none, "Introduction")

We should start with Chapters 1-5, then Chapter 8, then Chapter 11. Try all the exercises in the
book.

= The Integers

#problem[Find an example illustrating why the hypothesis that $a != 0$ is necessary in the statement of Lemma 1.2.]
#solution[We have a trivial example. If $a = 0$ and $b = 1$, then $b divides a$, as $a = 0b$.]

#problem[Let $n$ be a positive integer, and let $a$ and $b$ be integers. Prove that $a$ and $b$ have the same remainder when divided by $n$ if and only if $a - b = n k$ for some integer $k$.]

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
#problem[Prove that if $R$ is an integral domain, then $(x)$ is a prime ideal in $R[x]$.]
#solution[By the first isomorphism theorem, $R[x]\/(x) tilde.equiv R$, $R[x]\/(x)$ is an integral domain, and $(x)$ is a prime ideal of $R[x]$.]
#problem[]
#problem[]
#problem[Prove that the ideal $(y+1,x^2+1)$ is maximal in $RR[x,y]$ and contains the ideal $(y-x^2)$. This ideal does not correspond to a point with real coordinates on the parabola $y-x^2$. (This phenomenon is one reason why algebraic geometry is ‘easier’ over the complex numbers.)]

#counter(heading).update(7)

= Modules and Abelian Groups

- A vector space over a field $k$ is a set $V$ endowed with $(+, dot)$ by scalars
- The polynomial ring $k[x]$ over $k$ is a $k$-vector space
- This includes the ring addition axioms as well as:
  - $forall v in V, 1_k v = v$
  - $forall a, b in k, forall v in V, (a b)v = a (b V)$
  - $ forall a, b in k, forall v, w in V, (a + b)v  &= a v + b v \ a( v + w) &= a v + a w $
- Scalar multiplication defines an _action_ of $k$ on $V$ (since these are not really the ring axioms, but they do look very similar)
- We see ideal absorption $(forall r in R)(forall a in I), a r in I$ is similar to scalar multiplication:\ $R times I arrow.r I$
- Think of ideal over ring as some vector space
  - We define an $R$-module as this generalization, notated $M$
- A vector space is a "$k$-module" (my words)—an $R$-module where $R$ is a field $k$
- Rings are (obviously) $R$-modules
- Define a ring homomorphism $f: R arrow.r S$. Define an action of $R$ on $S$ by  $forall r in R, forall s in S, r dot s = f(r)s$.
  - This is an $R$-module over $S$
- $R\/I$ is also an $R$-module. $r dot (1 + I)$ generates $R\/I$
  - Then this is a cyclic $R$-module
- If ground ring is $ZZ$ then we can have the action $n dot a = underbrace((1+ dots.c + 1), n "times") dot a = underbrace(1 dot a + dots.c + 1 dot a, n "times")$


#counter(heading).update(10)

= Groups—Preliminaries