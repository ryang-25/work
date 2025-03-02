#import "@preview/great-theorems:0.1.2": *
#import "@preview/headcount:0.1.0": *

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
  )
  set text(size: 12pt)
  set heading(numbering: "1.")

  // Set common enum numbering and add some spacing.
  set enum(numbering: "a)", indent: 1em)
  set list(indent: 1em)

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

#let tcounter = counter("counter")
#let problem = mathblock(
  blocktitle: "Problem",
  counter: tcounter,
)
#let solution = proofblock(
  blocktitle: "Solution",
  prefix: [_Solution._],
  suffix: []
)

#let lifblock(body) = {
  set par(hanging-indent: 2.5em)
  h(0.5em)
  $==>$
  h(0.5em)
  body
}

#let rifblock(body) = {
  set par(hanging-indent: 2.5em)
  h(0.5em)
  $<==$
  h(0.5em)
  body
}
