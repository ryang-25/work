#!/usr/bin/env sage
"""
Some code from some exercises
"""

from itertools import combinations
from functools import reduce
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

def snd_dev_test(slns, fxx, fyy, fxy):
  out = []
  for x, y in slns:
    fxx_ = fxx(x, y)
    d = (fxx_*fyy - fxy^2)(x=x, y=y)
    if d > 0:
      if fxx_ > 0:
        out.append("min")
      elif fxx_ < 0:
        out.append("max")
      else:
        out.append("indeterminate")
    elif d < 0:
      out.append("saddle")
    else:
      raise ValueError
  return out

# Exercise 5
def s28_e5():
  print("Exercise 5")
  P.<x,y> = QQ[]
  f = (x^2 + y^2 - 4)*(x^2 + y^2 - 1) + (x - 3/2)^2 + (y - 3/2)^2
  
  fx = f.derivative(x)
  fxx = f.derivative(x, 2)
  fy = f.derivative(y)
  fyy = f.derivative(y, 2)

  fxy = fx.derivative(y)

  I = ideal(fx, fy)
  B = I.groebner_basis()
  B = [SR(g) for g in B] # convert to symbolic ring
  solutions = solve(B, var("x y"))
  solutions = format_solutions(solutions, SR)
  maxima = snd_dev_test(solutions, fxx, fyy, fxy)
  pretty_print(solutions)
  print(maxima)

# Exercise 10
def s28_e10():
  print("Exercise 10")
  P.<λ,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x - 1 - 2*λ*x^3, y - 1 - λ*y, z - 1 - λ*z, x^4 + y^2 + z^2 - 1)
  E = I.groebner_basis()[-1].univariate_polynomial()
  root = find_root(E, a=0, b=1) # find the positive root
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

# Section 3.1

# Exercise 2
def s31_e2():
  print("Exercise 2")
  def run():
    I = ideal(x^2 + 2*y^2 - 3, x^2 + x*y + y^2 - 3)
    # eliminate
    E = I.groebner_basis()[-1]
    v_type = E.variable(0)
    print("polynomial", E)
    rts = E.univariate_polynomial().roots(multiplicities=False)
    print("roots", rts)
    stms = [(I.basis.subs({v_type: v}), v) for v in rts]
    for stm, v in stms:
      stm = [SR(s) for s in stm]
      is_y = str(E.variable(0)) == "y"
      sol_v = var("x") if is_y else var("y")
      sls = [eq[0].rhs() for eq in solve(stm, sol_v)]
      out = (sls, v) if str(E.variable(0)) == "y" else (v, sls)
      print("solutions", *out, sep=", ")

  P.<x,y> = PolynomialRing(QQ, order="lex")
  run()
  P.<y,x> = PolynomialRing(QQ, order="lex")
  run()
  print()

def poly_roots(ps):
  res = []
  for p in ps:
    res.extend(p.univariate_polynomial().roots(multiplicities=False))
  out = set()
  for i in range(len(res)):
    for j in range(i, len(res)):
      if math.isclose(res[i], res[j], rel_tol=1e-15):
        out.add(round(res[i], 15))
  return out

def s31_e3():
  print("Exercise 3")
  def run(R):
    I = ideal(x^2 + 2*y^2 - 2, x^2 + x*y + y^2 - 2)
    E = I.groebner_basis()[-1]
    E_U, T = E.univariate_polynomial(), E.variable(0)
    roots = E_U.roots(ring=CC, multiplicities=False)
    # change to complex
    R = R.change_ring(CC)
    B = I.change_ring(R).basis

    solutions = []
    for root in roots:
      C = B.subs({T: root})
      t = C.variables()[0]
      solutions.extend([{t: p, T: root} for p in poly_roots(C)])
    print(f"Roots of {E}:", roots)
    print(f"Rational roots:", E_U.roots(ring=QQ, multiplicities=False))
    print("Solutions:", solutions)
  
  P.<x,y> = PolynomialRing(QQ, order="lex") # start with rationals
  run(P)


def s31_e4():
  print("Exercise 4")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x^2 + y^2 + x^2 - 4, x^2 + 2*y^2 - 5, x*z - 1)
  print("For the ideal", I.basis)
  B = I.groebner_basis()
  y_g = B[-2].univariate_polynomial()
  assert y_g.roots() == []
  print("No rational roots found for", y_g)
  print()


def s31_e7():
  print("Exercise 7")
  P.<t,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(t^2 + x^2 + y^2 + z^2, t^2 + 2*x^2 - x*y - z^2, t + y^3 - z^3)
  B = I.groebner_basis()
  Q.<x,y,z> = PolynomialRing(QQ, order="degrevlex")
  G = [Q(f) for f in B if f in Q]
  J = ideal(G)
  C = J.groebner_basis()
  T = TermOrder("elim", force=True)
  R.<t,x,y,z> = PolynomialRing(QQ, order=T)
  I = ideal(t^2 + x^2 + y^2 + z^2, t^2 + 2*x^2 - x*y - z^2, t + y^3 - z^3)
  K = [R(f) for f in C] + [t + y^3 - z^3]

  print("First, generators of the ideal with lex order", G)
  print()
  print("Then the reduced generators under degrevlex are", C, sep="\n")
  print()
  print("If we add t + y^3 - z^3 to the set", K, sep="\n")
  print("And the generators under elim order are", I.groebner_basis())

def s32_e2():
  print("Exercise 2")
  P.<x,y> = PolynomialRing(CC, order="lex")
  I = ideal(y*x^3+x^2, y^3*x^2+y^2, y*x^4 + x^2 + y^2)
  B = I.groebner_basis()
  print("The Groebner basis of I is", B)
  assert B[-1] == y^2 # a
  S = Sequence(y, y^3, y)
  I_ = ideal(I.basis + S)
  B_ = I_.groebner_basis()
  print("The Groebner basis of I' is", B_)
  G = []
  for g in I.basis:
    G.extend([g-g.lt(), g.lt()])
  I__ = ideal(G)
  B__ = I__.groebner_basis()
  print("The Groebner basis of tilde I is", B__)
  print()

def s32_e3():
  print("Exercise 3")
  P.<x,y,z> = PolynomialRing(CC, order="lex")
  I = ideal(x^2 + y^2+ z^2 + 2, 3*x^2 + 4*y^2 + 4*z^2 + 5)
  B = I.groebner_basis()
  print("The Groebner basis of I is", B)

def s33_e6():
  print("Exercise 6")
  P.<u,v,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x - u*v, y - u^2, z - v^2)
  B = I.groebner_basis()
  Q.<v,x,y,z> = QQ[]
  R.<x,y,z> = QQ[]
  V2_B = [g for g in B if g in R]
  I1 = [g for g in B if g in Q]
  assert I1[0].lt() == v^2
  assert B[0].lt() == u^2
  print("The Groebner basis of I is", B)
  print("The equation of the variety is", V2_B)
  print("The first elimination ideal is", I1)
  print("Then the leading coefficient is 1, so it extends to v")
  print("Then the leading coefficient is 1, so it extends to u")
  print()

def s33_e7():
  print("Exercise 7")
  P.<u,v,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x - u*v, y - u*v^2, z - u^2)
  B = I.groebner_basis()
  Q.<v,x,y,z> = QQ[]
  R.<x,y,z> = QQ[]
  V2_B = [g for g in B if g in R]
  I1 = [g for g in B if g in Q]
  print("The Groebner basis of I is", B)
  print("The equation of the variety is", V2_B)
  print("The first elimination ideal is", I1)
  print()

def find_corollary_4(seq, var):
  for p in seq:
    if p.lt() == 0:
      pass

def s33_e8():
  print("Exercise 8")
  P.<u,v,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x - (3*u + 3*u*v^2 - u^3), y - (3*v + 3*u^2*v - v^3), z - (3*u^2 - 3*v^2))
  B = I.groebner_basis()
  Q.<v,x,y,z> = QQ[]
  R.<x,y,z> = QQ[]
  V2_B = [g for g in B if g in R]
  I1 = [g for g in B if g in Q]
  assert v^3 + 1/2*v*z + 3/2*v - 1/2*y in I1
  assert u^2 - v^2 - 1/3*z in I
  # print("The Groebner basis of I is", *B, sep="\n\n")
  print("The equation of the variety is", V2_B)
  # print("The first elimination ideal is", *I1, sep="\n\n")
  print("Then the leading coefficient is 1, so it extends to v")
  print("Then the leading coefficient is 1, so it extends to u")
  print()

def s33_e9():
  print("Exercise 9")
  P.<u,v,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x - u*v, y - v, z - u^2)
  B = I.groebner_basis()
  Q.<v,x,y,z> = QQ[]
  R.<x,y,z> = QQ[]
  V2_B = [g for g in B if g in R]
  I1 = [g for g in B if g in Q]
  assert I1[0].lt() == v
  assert B[0].lt() == u^2
  print("The Groebner basis of I is", B)
  print("The equation of the variety is", V2_B)
  print("The first elimination ideal is", I1)
  print("Extends over the complex numbers.")

def s33_e14():
  print("Exercise 14")
  P.<t,x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x*(1 + t^3) - 3*t,y*(1 + t^3) - 3*t^2)
  B = I.groebner_basis()
  Q.<x,y,z> = QQ[]
  I1 = [g for g in B if g in Q]

  print("The Groebner basis of I is", B)
  # print("The equation of the variety is", V2_B)
  print("The first elimination ideal is", I1)
  # print("Extends over the complex numbers.")
  print()


def s35_e4():
  print("Exercise 4")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(x^2*y + x*z + 1, x*y - x*z^2 + z - 1)
  B = I.groebner_basis()
  print(B)

def s41_e1():
  print("Exercise 1")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  I = ideal(y - x^2, z - x^3)
  J = ideal((y - x^2)^2 + (z - x^3)^2)
  print(I.groebner_basis())
  print(J.groebner_basis())
  print()

def s42_e11():
  print("Exercise 11")
  P.<x,y> = QQ[] # 1 poly field doesn't support radical calculation.
  f1 = x^5 - 2*x^4 + 2*x^2 - x
  f2 = x^5 - x^4 - 2*x^3 + 2*x^2 + x - 1
  I = ideal([f1, f2])
  # print(I.variety())
  print(I.radical().basis)

def s43_e8():
  print("Exercise 8")
  P.<x,y,z> = PolynomialRing(QQ, order="lex")
  f = x^4 + x^3*y + x^3*z^2 - x^2*y^2 + x^2*y*z^2 - x*y^3 - x*y^2*z^2 - y^3*z^2
  I = ideal(f)
  g = x^4 + 2*x^3*z^2 - x^2*y^2 + x^2*z^4 - 2*x*y^2*z^2 - y^2*z^4
  J = ideal(g)
  K = I.intersection(J)
  R = (I*J).radical()
  gcd_ = f.gcd(g)

  p = x^2 + x*y + x*z + y*z
  q = x^2 - x*y - x*z + y*z
  I = ideal(f, g)
  J = ideal(p, q)
  L = I.intersection(J)

  print("I intersect J", K.basis)
  print("Radical of IJ", R.basis)
  print("gcd(f,g)", gcd_)
  print("(f, g) intersect (p, q)", L.basis)


def s46_e5():
  print("Exercise 5")
  P.<z,y,x> = PolynomialRing(QQ, order="lex")
  J = ideal([x*z - y^2, x^3 - y*z, z^2 - x^2*y])
  B = J.groebner_basis()
  print(B)

def s51_e2():
  print("Exercise 2")
  P.<t,u,v,x,y,z> = PolynomialRing(QQ, order="lex")
 #I = ideal([x-t^2, y-t, z-t^3, x^2-y-(t^2-t), y^2-t^2])

  # I = ideal([x-t^2, y-t, z-t^3, x*y-u, z+x^2*y^2-v])

  I = ideal([x-t^2, y-t, z-t^3, x*y-u, z+x^2*y^2-v])
  B = I.groebner_basis()
  print(B)

def s53_e5():
  P.<x,y> = PolynomialRing(QQ, order="lex")
  B = Sequence([x^2 + y - 1, x*y - 2*y^2 + 2*y, y^3 - 7/4*y^2 + 3/4*y])
  l = [1, x, y, y^2]
  for f, g in zip(l, l):
    if f*g == 1:
      continue
    print(f"{f} * {g}: {multi_rem(f*g, B)}")
  for f, g in combinations(l, 2):
    print(f"{f} * {g}: {multi_rem(f*g, B)}")
  print()


def s53_e6():
    P.<x5,x4,x3,x2,x1> = PolynomialRing(QQ, order="lex")
    I = ideal([x3-x1^2, x4-x1*x2, x2*x4 - x1*x5, x4^2 - x3*x5])
    B = I.groebner_basis()
    print(B)

def main():
  """
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
  s23_e6()

  print("\nSection 2.7")
  s27_e2()

  print("\nSection 2.8")
  s28_e3()
  s28_e5()
  s28_e10()
  s28_e11()
  """
  """
  s31_e2()
#  s31_e3()
  s31_e4()
  s31_e7()
  print("\nSection 3.2")
  s32_e2()
  s32_e3()

  print("\nSection 3.3")
  s33_e6()
  s33_e7()
  s33_e8()
  s33_e9()
  s33_e14()


  s35_e4()

  s41_e1()

  print("\nSection 4.2")
  s42_e11()

  print("\nSection 4.3")
  s43_e8()

  print("\nSection 4.5")
  s46_e5()
  """

  print("\nSection 5.1")
  s51_e2()

  s53_e5()
  s53_e6()

if __name__ == '__main__':
  main()
