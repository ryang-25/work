#import "@preview/ctheorems:1.1.3": *
#show: thmrules

#let homework(
  doc,
  assign: none,
  date: datetime,
  name: none,
) = {
  set text(font: "New Computer Modern", size: 12pt)
  set page(
    margin: (top: 1.5in, left: 1in, right: 1in, bottom: 1in),
    header: [
      #set text(weight: "bold")
      #set par(leading: 0.75em)
      CSC 403
      #h(1fr) Name: #text(style: "italic", weight: "regular")[#name] \
      #assign \
      #date.display("[day] [month repr:long] [year]")
      #line(length: 100%, stroke: 2pt)
    ],
    header-ascent: 10%
  )
  doc
}

#let problem = thmbox(
  "problem",
  "Problem",
  inset: 0em
)

#show : homework.with(
  assign: "Final Project Report",
  date: datetime.today(),
  name: "Roland Yang"
)



== IMO 2020

#problem("IMO 2020", number: "2")[If ]


