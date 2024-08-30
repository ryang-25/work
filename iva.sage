#!/usr/bin/env sage

# Exercise 8, Section 1.4
x = PolynomialRing(QQ, "x").gen()
e8a = (x^4 + x^2 + 1).gcd(x^4 - x^2 - 2*x - 1).gcd(x^3-1)
e8b = (x^3 + 2*x^2 - x - 2).gcd(x^3 - 2*x^3 - x + 2).gcd(x^3 - x^2 - 4*x + 4)
assert e8a == x^2 + x + 1
assert e8b == x - 1

# Exercise 9
e9 = (x^3 + x^2 - 4*x - 4).gcd(x^3 - x^2 - 4*x+4).gcd(x^3 - 2*x^2 - x + 2)
assert e9 == x - 2

# Exercise 1, Section 1.1
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

# Exercise 1, Section 2.3
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
