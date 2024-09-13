#!/usr/bin/env sage
"""
Some code from some exercises
"""

from itertools import combinations
from functools import reduce
from sage.quadratic_forms.qfsolve import qfsolve, qfparam
from sage.rings.polynomial.toy_buchberger import spol

# Section 1.4

def s14_e8():
  print("Exercise 8")
  P.<x> = QQ[]
  ra = reduce(gcd, [x^4 + x^2 + 1, x^4 - x^2 - 2*x - 1, x^3-1])
  rb = reduce(gcd, [x^3 + 2*x^2 - x - 2, x^3 - 2*x^3 - x + 2, x^3 - x^2 - 4*x + 4])
  assert ra == x^2 + x + 1
  assert rb == x - 1

def s14_e9():
  print("Exercise 9")
  P.<x> = QQ[]
  # Exercise 9
  r = reduce(gcd, [x^3 + x^2 - 4*x - 4, x^3 - x^2 - 4*x+4, x^3 - 2*x^2 - x + 2])
  assert r == x - 2

# Section 2.1

def s21_e1():
  print("Exercise 1")
  P.<x> = QQ[]
  c = (x^4 - 6*x^2 + 12*x - 8).gcd(2*x^3-10*x^2 + 16*x - 8)
  d = (x^9 - 1).gcd(x^5 + x^3 - x^2 - 1)
  assert c == 1
  assert d == x^3 - 1

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

# Section 2.3

def s23_e1_3():
  print("Exercise 1 & Exercise 3")
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

def s23_e5():
  print("Exercise 5")
  P.<x,y,z> = PolynomialRing(QQ, order="deglex")
  f = x^3 - x^2*y - x^2*z + x
  f1 = x^2*y - z
  f2 = x*y - 1
  r1 = multi_rem(f, [f1, f2])
  r2 = multi_rem(f, [f2, f1])
  assert r1 != r2
  assert multi_rem(r1 - r2, [f1, f2]) == x - z

def s23_e6():
  print("Exercise 6")
  # something here?

def spoly(f, g):
  m = f.lm().lcm(g.lm())
  return m//f.lt() * f - m//g.lt() * g

# Section 2.6

def s26_e6():
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
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

# Section 2.7

def s27_e2():
  print("Exercise 2")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  IA = ideal(x^2*y - 1, x*y^2 - x)
  BA = toy_groebner_basis(IA.basis)

  IB = ideal(x^2 + y, x^4 + 2*x^2*y + y^2 + 3)
  BB = toy_groebner_basis(IB.basis)

  IC = ideal(x - z^4, y - z^5)
  BC = toy_groebner_basis(IC.basis)

  P.<x,y,z> = PolynomialRing(QQ, order="deglex")
  IA = ideal(x^2*y - 1, x*y^2 - x)
  BA = toy_groebner_basis(IA.basis)

  IB = ideal(x^2 + y, x^4 + 2*x^2*y + y^2 + 3)
  BB = toy_groebner_basis(IB.basis)

  IC = ideal(x - z^4, y - z^5)
  BC = toy_groebner_basis(IC.basis)

# Section 2.8

def s28_e1():
  print("Exercise 1")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(-x^3 + y, x^2*y - z)
  f = x*y^3 - z^2 + y^5 - z^3
  assert f in I

def s28_e2():
  print("Exercise 2")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x*z - y, x*y + 2*z^2, y - z)
  f = x^3*z - 2*y^2
  assert not f in I

def format_solutions(slns, ring):
  """
  There's definitely a function for this.
  """
  return [tuple([ring(v.rhs()) for v in s]) for s in slns]

def s28_e3():
  print("Exercise 3")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x^2 + y^2 + z^2 - 1, x^2 + y^2 + z^2 - 2*x, 2*x - 3*y - z)
  B = I.groebner_basis()
  B = [SR(g) for g in B] # convert to symbolic ring
  solutions = solve(B, var("x y z"))
  pretty_print(format_solutions(solutions, SR))

# Exercise 5
def s28_e5():
  print("Exercise 5")
  P.<x,y> = QQ[]
  f = (x^2 + y^2 - 4)*(x^2 + y^2 - 1) + (x - 3/2)^2 + (y - 3/2)^2
  
  fx = f.derivative(x)
  fxx = fx.derivative(x)
  fy = f.derivative(y)
  fyy = fy.derivative(y)

  I = ideal(fx, fy)
  B = I.groebner_basis()
  B = [SR(g) for g in B] # convert to symbolic ring
  solutions = solve(B, var("x y"))
  pretty_print(format_solutions(solutions, SR))

# Exercise 10
def s28_e10():
  print("Exercise 10")
  P.<λ,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x - 1 - 2*λ*x^3, y - 1 - λ*y, z - 1 - λ*z, x^4 + y^2 + z^2 - 1)
  B = I.groebner_basis()
  P.<z> = QQ[]
  root = find_root(P(B[-1]), a=0, b=1) # find the positive root
  # Solution with hand-calculated Lagrange multipliers
  x = var("x")
  y = z = root
  x = find_root(x^4 + y^2 + z^2 - 1, a=0, b=1)
  print("Solution:", end=" ")
  print(x, y, z, sep=", ")

# Exercise 11
def s28_e11():
  print("Exercise 11")
  P.<a,b,c> = QQ[]
  I = ideal(a+b+c-3, a^2+b^2+c^2-5, a^3+b^3+c^3-7)
  f = lambda n : I.reduce(a^n + b^n + c^n)
  print(list(f(r) for r in range(4, 9)))


def main():
  print("Section 1.4")
  s14_e8()
  s14_e9()

  print("\nSection 2.1")
  s21_e1()

  print("\nSection 2.3")
  s23_e1_3()
  s23_e5()
  s23_e6()

  print("\nSection 2.6")

  print("\nSection 2.7")
  s27_e2()

  print("\nSection 2.8")
  s28_e3()
  s28_e5()
  s28_e10()
  s28_e11()

if __name__ == '__main__':
  main()