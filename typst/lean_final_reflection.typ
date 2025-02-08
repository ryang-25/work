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

#let lemma = thmbox(
  "lemma",
  "Lemma",
  titlefmt: emph,
  inset: 0em).with(numbering: none)
#let problem = thmbox(
  "problem",
  "Problem",
  inset: 0em
)
#let solution = thmproof("solution","Solution",inset: 0em)
#show: thmrules.with(qed-symbol: "")


#show: homework.with(
  assign: "Final Project ",
  date: datetime.today(),
  name: "Roland Yang"
)

This is my own final project report.

To choose a proof to tackle, I went through the AoPS site and looked for "simple algebra or number theory problems. I
then went and found solutions published by IMO themselves to guide my formalization. I ended up producing partial
proofs for both of my problems.

== IMO 2024

#problem("IMO 2024", number: "2")[Determine all pairs of numbers $(a, b)$ of positive integers for which there exist
 integers $g$ and $N$ such that $ gcd(a^n + b, b^n + a) = g $ holds for all integers $n >= N$.]

#solution[We show that the only solution is $(a,b) = (1,1)$.

We can verify that it is true when $g = 2$
#lemma[We have that $g = gcd(a,b)$ or $g = 2gcd(a,b)$]]

== IMO 2022

_Solution:_ #link("https://www.imo-official.org/problems/IMO2022SL.pdf")[www.imo-official.org/problems/IMO2022SL.pdf]

#problem("IMO 2022", number: "5")[Find all triples of positive integers $(a,b,p)$ with $p$ prime and $ a^p = b! + p. $]

#solution[$(2,2,2)$ and $(3,4,3)$.

We can break the proof into inequalities. We can see that $a > 1$ trivially,
with

```lean
have gt_one : 1 < a := by
  apply lt_of_le_of_ne ha
  by_contra! hc
  replace hc : 1 = a := hc
   [← hc, one_pow] at h
  have := factorial_pos b
  have := hp.pos
  have : 1 < a := by omega
  exact this.ne.elim hc
```

Then we split on three cases: when $a < p$, when $a > p$, or when $a = p$.

#grid(
  columns: (1fr, auto),
  row-gutter: 0.75em,
  [- $a < p$:], [We have either $a <= b$ or $a > b$ \
  For $a <= b arrow.r a! <= b! $ so $a | a^p - b!$ and $a | p$, but $a eq.not 1 eq.not p$, a contradiction.],
  [- $a > p$],[],
  [- $a = p$], [],
)]
