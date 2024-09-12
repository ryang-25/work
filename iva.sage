#!/usr/bin/env sage
"""
I really wish this was better
"""

from itertools import combinations
from functools import reduce
from sage.quadratic_forms.qfsolve import qfsolve, qfparam
from sage.rings.polynomial.toy_buchberger import spol

# Exercise 8, Section 1.4
P.<x> = QQ[]
e8a = reduce(gcd, [x^4 + x^2 + 1, x^4 - x^2 - 2*x - 1, x^3-1])
e8b = reduce(gcd, [x^3 + 2*x^2 - x - 2, x^3 - 2*x^3 - x + 2, x^3 - x^2 - 4*x + 4])
assert e8a == x^2 + x + 1
assert e8b == x - 1

# Exercise 9
e9 = reduce(gcd, [x^3 + x^2 - 4*x - 4, x^3 - x^2 - 4*x+4, x^3 - 2*x^2 - x + 2])
assert e9 == x - 2

# Exercise 1, Section 2.1
e1c = (x^4 - 6*x^2 + 12*x - 8).gcd(2*x^3-10*x^2 + 16*x - 8)
assert e1c == 1

e1d = (x^9 - 1).gcd(x^5 + x^3 - x^2 - 1)
assert e1d == x^3 - 1

def multi_rem(p, ds):
  rem = 0
  while p != 0:
    divisible = False
    for d in ds:
      if p.lt()//d.lt():
        p -= (p.lt()//d.lt())*d
        divisible = True
        break
    if not divisible:
      rem += p.lt()
      p -= p.lt()
  return rem

# Exercises 1, 3, Section 2.3
P.<x,y> = PolynomialRing(QQ, order="deglex")
f = x^7 * y^2 + x^3*y^2 - y + 1
f1 = x*y^2 - x
f2 = x - y^3
assert multi_rem(f, [f1, f2]) == multi_rem(f, [f2, f1]) == x^7 + x^3 - y + 1

P.<x,y> = PolynomialRing(QQ, order="lex")
f = x^7*y^2 + x^3*y^2 - y + 1
f1 = x*y^2 - x
f2 = x - y^3
assert multi_rem(f, [f1, f2]) == 2*y^3 - y + 1
assert multi_rem(f, [f2, f1]) == y^23 + y^11 -y + 1

# Exercise 5
P.<x,y,z> = PolynomialRing(QQ, order="deglex")
f = x^3 - x^2*y - x^2*z + x
f1 = x^2*y - z
f2 = x*y - 1
r1 = multi_rem(f, [f1, f2])
r2 = multi_rem(f, [f2, f1])
assert r1 != r2
assert multi_rem(r1 - r2, [f1, f2]) == x - z

# Exercise 6
# something here?

def spoly(f, g):
  m = f.lm().lcm(g.lm())
  return m//f.lt() * f - m//g.lt() * g

# Exercise 6, Section 2.6
P.<x,y,z> = PolynomialRing(QQ, order="lex")
# Exercise 6, Section 2.6
f, g = 4*x^2*z - 7*y^2, x*y*z^2 + 3*x*z^4
assert spol(f,g) == -3*x^2*z^4 - 7/4*y^3*z # a
f, g = x^4*y - z^2, 3*x*z^2 - y
assert spol(f,g) == x^3*y^2/3 - z^4 # b
f, g = x^7*y^2*z + 2*I*x*y*z, 2*x^7*y^2*z + 4
assert spol(f,g) == 2*I*x*y*z - 2 # c
f, g = x*y + z^3, z^2 - 3*z
assert spol(f,g) == 3*x*y*z + z^5 # d

def toy_groebner_basis(fs: Sequence):
  """
  interactive groebner basis calculation
  """
  g = set(fs)
  g_ = {}
  while g != g_:
    g_ = g.copy()
    for p, q in combinations(g_, 2):
      s = spol(p, q)
      r = multi_rem(s, g_)
      if r != 0:
        g.add(r)
  return Sequence(g)

# Exercise 2, Section 2.7
P.<x,y,z> = PolynomialRing(QQ, order="lex")
ia = ideal(x^2*y - 1, x*y^2 - x)
ba = toy_groebner_basis(ia.basis)
print(ba)
assert set(ba) == {x^2*y - 1, y^3 - y, x^2 - y, y^2 - 1, x*y^2 - x}
ib = ideal(x^2 + y, x^4 + 2*x^2*y + y^2 + 3)
bb = set(toy_groebner_basis(ib.basis))
# assert bb == {}
ic = ideal(x - z^4, y - z^5)
bc = toy_groebner_basis(ic.basis)
assert set(bc) == {x - z^4, y - z^5}

P.<x,y,z> = PolynomialRing(QQ, order="deglex")
ia = ideal(x^2*y - 1, x*y^2 - x)
ba = toy_groebner_basis(ia.basis)
assert set(ba) == {x^2*y - 1, y^3 - y, x^2 - y, y^2 - 1, x*y^2 - x}
ib = ideal(x^2 + y, x^4 + 2*x^2*y + y^2 + 3)
bb = toy_groebner_basis(ib.basis)
print(bb)
ic = ideal(x - z^4, y - z^5)
bc = toy_groebner_basis(ic.basis)
print(bc)

# Exercise 1, Section 2.8
P.<x,y,z> = PolynomialRing(QQ, order="lex")
i = ideal(-x^3 + y, x^2*y - z)
f = x*y^3 - z^2 + y^5 - z^3
assert f in i

# Exercise 2
i = ideal(x*z - y, x*y + 2*z^2, y - z)
f = x^3*z - 2*y^2
assert not f in i

def format_solutions(slns):
  """
  There's definitely a function for this.
  """
  return [tuple([v.rhs() for v in s]) for s in slns]

# Exercise 3
i = ideal(x^2 + y^2 + z^2 - 1, x^2 + y^2 + z^2 - 2*x, 2*x - 3*y - z)
b = i.groebner_basis()
b = [SR(g) for g in b]
# solve(b, var("x"), var("y"), var("z"))

# Exercise 5
f = (x^2 + y^2 - 4)*(x^2 + y^2 - 1) + (x - 3/2)^2 + (y - 3/2)^2
fx = f.derivative(x)
fy = f.derivative(y)
i = ideal(fx, fy)
b = i.groebner_basis()
b = [SR(g) for g in b]
q = solve(b, var("x"), var("y"))
print(format_solutions(q))

# Exercise 10

# Exercise 11
