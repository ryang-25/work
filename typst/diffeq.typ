#import "notes.typ": *
#show: great-theorems-init

// template
#let notes(
  title: none,
  subject: none,
  author: [],
  date: datetime.today(),
  doc
) = {
  set page(
    paper: "us-letter",
    header: [
      #set block(spacing: 0.5em)
      #title #h(1fr) #text(style: "italic", subject)
      #line(length: 100%)
    ],
    footer: context [
     #align(center, counter(page).display("1 / 1", both:true))
    ],
    margin: 1in
  )
  set text(size: 12pt)
  // set heading(numbering: "1.")
  // show heading.where(level: 3): set heading(numbering: none)

  // front matter
  align(right,
    block(
      stroke: (left: 1pt),
      inset: (left: 0.75em, right: 1.5em, top: 1.25em, bottom: 2em)
    )[
    #set align(left)
    #set text(size: 13pt)
    #text(size: 28pt, weight: "bold", title)
    #v(-0.5em)
    #date.display("[day] [month repr:long] [year]") \
    #author
  ])
  doc
}

#show: notes.with(
  title: "Differential Equations",
  subject: "Class Notes",
  author: "Roland Yang",
)

#let dy = $dif y$
#let dx = $dif x$

// #theorem[Hi]

#heading(numbering: none, "Introduction")

= #datetime(month: 2, day: 7, year: 2025).display("[month]/[day]/[year repr:last_two]")

#let du = $dif u$

We're talking about change of variables / substitution today. We already know how to solve first order linear
differential equations of the form $ dy/dx = a y + b $ where $a, b$ are functions of $x$. But we want to use
substitution to solve other equations. First, we recognize that we can either substitute the dependent variable $y$ or
the independent variable $x$. Take the below problem as example:

#problem[Solve the differential equation $ dy/dx = y/2 + x/y $]

We don't have the tools to solve something like this yet, but let's substitute $y$. Introduce the variable $u = y^2$ and solve.

Why do we use $u = y^2$? If we look at a power $n$ where $y = u^n$, then using the notation of partial differentiation$ y_x &= n u^(n-1) u_x \ n u_x &= a u + b u^(1-n) $ then our choice of $n$ depends on the power of $b u^(1-n)$. Going back to the original problem... TBD

For the below problem, we have an isocline $y_x|_(y = x) = 0$ except obviously at the singular point

#problem[$(y^2 - x^2)dx + x y dy = 0$]