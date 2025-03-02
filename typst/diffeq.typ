#import "notes.typ": *
#show: great-theorems-init

#show: notes.with(
  title: "Differential Equations",
  subject: "Class Notes",
  author: "Roland Yang",
)
#show heading: reset-counter(tcounter, levels: 2)

#let date_display(month, day, year) = datetime(month: month, day: day, year:
  year).display("[month]/[day]/[year repr:last_two]")

#let dy = $dif y$
#let dx = $dif x$
#let dz = $dif z$

#heading(numbering: none, "Introduction")

= #date_display(2, 7, 2025)


#let du = $dif u$

== Bernoulli Equations

We're talking about change of variables / substitution today. We already know
how to solve first order linear differential equations of the form $ dy/dx =
  a y + b $ where $a, b$ are functions of $x$. But we want to use substitution
to solve other equations. First, we recognize that we can either substitute the
dependent variable $y$ or the independent variable $x$. Take the below problem
as example:

#problem[Solve the differential equation $ dy/dx = y/2 + x/y $]

We don't have the tools to solve something like this yet, but let's substitute $y$. Introduce the variable $u = y^2$ and solve.

Why do we use $u = y^2$? If we look at a power $n$ where $y = u^n$, then using
the notation of partial differentiation$ y_x &= n u^(n-1) u_x \ n u_x &=
  a u + b u^(1-n) $ then our choice of $n$ depends on the power of $b u^(1-n)$.
Going back to the original problem... TBD

For the below problem, we have an isocline $y_x|_(y = x) = 0$ except obviously at the singular point

#problem[$(y^2 - x^2)dx + x y dy = 0$]

= #date_display(2, 11, 2025)

Today we're covering _exact equations_ (2.4), which are of the form
$M dx + N dy = 0$. Take $ 3x^2 - x y + y^2 = 0 $ and call it $f$. From
multivariable, we know the differential $dz = f_x dx + f_y dy$ where $z = f(x,y)$
so taking the partials $ f_x &= 6x - y \ f_y &= 2y - x $ we can construct a
differential equation $(6x - y)dx + (2y - x)dy = 0$ where $f_(x y) = f_(y x) = -1$
so it is in _exact form_. Stating this more formally: if $ M_y = N_x $
*it is in exact form*.

#problem[Check the exactness of $\(y^2 - y/(2sqrt(x))\)dx + \(2 x y - sqrt(x) + 1\)dy$.]
#solution[$ M_y = N_x = 2y - 1/(2sqrt(x)) $ We can also go further, and see that
$ f = integral f_x dx = x y^2 - y sqrt(x) + C(y) $ and since $N_y = 2 x y -
  sqrt(x) + 1$ we see that $C_y = 1$ so $C$ must be a linear function in
$C = y + C_1$.

We can then go further and choose points of relation $f = x y^2 -
  y sqrt(x) + y + C$ and draw the level curves as well.
]

We then move on to using the integrating factor $mu$ to get non-exact equations
in the exact form.
#tcounter.update(45)
#problem[Use the integrating factor to solve $y(2e^x + 4x)dx + 3(e^x + x^2)dy = 0$]
#solution[
  Quickly verify that $M_y eq.not N_x$, where $M_y = 2e^x + 4x$ and
  $N_x = 3\(e^x + 2x\)$. Then we try multiplying as in 2.2 by a integrating factor
  $mu$. We have three cases: where $mu = mu(x)$, $mu = mu(y)$, or $mu = mu(x,y)$.
  We first try $mu = mu(x)$, and obtain $ M_y &= mu\(2e^x + 4x\) \
    N_x &= 3\(e^x + 2x\) mu + 3\(e^x + x^2\) mu_x $
We can ask ourselves: is there a $mu$ that satisfies the equation
$
  mu\(2e^x + 4x\) &= 3\(e^x + 2x\) mu + 3\(e^x + x^2\) mu_x \
  mu\(2e^x + 4x - 3e^x - 6x\) &= 3\(e^x + x^2\) mu_x \
  integral (dif mu)/mu &= -1/3 integral (e^x + 2x)/(e^x + x^2) dx \
  ln abs(mu) &= -1/3 ln abs(e^x + x^2) + C \
  mu &= 1/root(3, e&x + x^2)
$

We can ignore the $C$ since we are searching for _a_ $mu$, not any.

]
