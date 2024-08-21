#!/usr/bin/env sage

# Exercise 8, Section 1.4
x = PolynomialRing(RationalField(), 'x').gen()
e8a = (x^4 + x^2 + 1).gcd(x^4 - x^2 - 2*x - 1).gcd(x^3-1)
e8b = (x^3 + 2*x^2 - x - 2).gcd(x^3 - 2*x^3 - x + 2).gcd(x^3 - x^2 - 4*x + 4)
print(e8a) # x^2 + x + 1
print(e8b) # x - 1

# Exercise 9
e9 = (x^3 + x^2 - 4*x - 4).gcd(x^3 - x^2 - 4*x+4).gcd(x^3 - 2*x^2 - x + 2)
print(e9) # x - 2
