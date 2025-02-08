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
  inset: 0em,
  padding: (top: 0em, bottom: -2em),
).with(numbering: none)
#let problem = thmbox(
  "problem",
  "Problem",
  inset: 0em
)
#let solution = thmproof(
  "solution",
  "Solution",
  inset: 0em,
  padding: (top: 0em, bottom: 0em),
)
#show: thmrules.with(qed-symbol: "")


#show: homework.with(
  assign: "Final Project Paper",
  date: datetime.today(),
  name: "Roland Yang"
)


== IMO 2022

_Solution:_ #link("https://www.imo-official.org/problems/IMO2022SL.pdf")[www.imo-official.org/problems/IMO2022SL.pdf]

#problem("IMO 2022", number: "5")[Find all triples of positive integers $(a,b,p)$ with $p$ prime and $ a^p = b! + p. $]

#solution[$(2,2,2)$ and $(3,4,3)$.

Clearly, $a > 1$.

#grid(
  columns: (1.75cm, auto),
  row-gutter: 0.75em,
  [*Case 1:*], [We have $a < p$. Then we either have $a <= b$ or $a > b$ which implies $a | a^p - b! = p$ leading to a contradiction, or $a < b$ which is also impossible since in this case we have $b! <= a! < a^p-p$, where the last inequality is true for any $p > a > 1$.],
  [*Case 2:*], [We have $a > p$. In this case $b! = a^p -p > p^p - p >= p!$ so $b > p$ which means that $a^p = b! + p$ is divisible by $p$. Hence $a$ is divisible by $p$ and $b! = a^p - p$ is not divisible by $p^2$. This means that $b < 2p$. If $a < p^2$ then $a/p < p$ divides both $a^p$ and $b!$ and hence it also divides $p = a^p - b!$ which is impossible. On the other hand, the case $a >= p^2$ is also impossible since then $a^p >= (p^2)^p > (2p - 1)! + p >= b! + p$.],
  [*Case 3:*], [We have $a = p$. In this case $b! = p^p - p$. Checking primes $p < 7$ gives either solutions or contradiction, so have $p >= 7$. Then by lifting the exponent $ v_2((p+1)!) <= v_2(b!) = v_2(p^(p-1)-1) =^(L T E) = v_2((p - 1)/2 dot (p + 1) dot (p - 1))$. On the RHS we have three factors of $(p + 1)!$, but due to $p + 1 >= 8$, there are at least 4 even numbers among $1, 2, dots, p + 1$ that divide $p + 1$, so this case is not possible.],
)]


== Code

```lean
-- todo: make this somewhat faster.

-- https://github.com/leanprover-community/mathlib4/blob/6c225037eb50056d12f06584a2dd4ee0008910d9/Archive/Imo/Imo1994Q1.lean#L16
-- https://www.imo-official.org/problems/IMO2022SL.pdf N4.

import Mathlib.Data.Nat.Factorial.Basic
import Mathlib.Data.Nat.Factorial.BigOperators
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.NumberTheory.Multiplicity
import Mathlib.Tactic.GCongr
import Mathlib.Tactic.IntervalCases
-- https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/Multiplicity.html

open Nat

def solutionSet : Set (ℕ × ℕ × ℕ) := {(2,2,2), (3,4,3)}

-- helpful!
lemma factorial_lt_pow {n : ℕ} : n ! ≤ n^n := by
  induction' n with n hn
  apply_rfl
  rewrite [← mul_factorial_pred n.add_one_pos, Nat.add_one_sub_one]
  calc (n + 1) * n !
      ≤ (n + 1) * n^n := by rel [hn]
    _ ≤ (n + 1) * (n + 1)^n := by rel [Nat.pow_le_pow_left n.lt_add_one.le _]
    _ = (n + 1)^(n + 1) := symm pow_succ'

-- a big thank you to Daniel Weber on Zulip for almost entirely writing this proof
-- omega works overtime
lemma factorial_lt_pow_sub {n : ℕ} (h: 2 < n): n ! < n^n - n := by
  have sub_one_pos : 0 < n - 1 := by omega
  have : (n - 2)^(n - 2) < n^(n - 2) := Nat.pow_lt_pow_left (by omega) (by omega)
  calc
    _ = n * (n - 1) * (n - 1 - 1)! := by rw [← mul_factorial_pred $ zero_lt_of_lt h,
      ← mul_factorial_pred sub_one_pos, mul_assoc]
    _ = n * (n - 1) * (n - 2)! := by apply_rfl
    _ ≤ n * (n - 1) * (n - 2)^(n - 2) := by rel [factorial_lt_pow]
    _ < n * (n - 1) * n^(n - 2) := by rel [this]
    _ = n^n - n^(n - 2) * n := by rw [Nat.mul_sub_left_distrib, mul_one,
      ← sq, mul_comm, Nat.mul_sub_left_distrib, Nat.pow_sub_mul_pow _ h.le]
    _ = n^n - n^(n - 1 - 1) * n := by apply_rfl
    _ = n^n - n^(n - 1 - 1) * n^1 := by rw [pow_one]
    _ = n^n - n^(n - 1) := by rw [Nat.pow_sub_mul_pow _ sub_one_pos.nat_succ_le]
    _ ≤ n^n - n := Nat.sub_le_sub_left (Nat.le_self_pow sub_one_pos.ne' _) _

lemma factorial_le_pow_sub {n : ℕ} (h: 2 ≤ n) : n ! ≤ n^n - n := by
  rcases h.eq_or_lt with rfl | hg
  apply_rfl
  exact (factorial_lt_pow_sub hg).le

-- needed for proof
lemma emultiplicity_div {p a b: ℕ} (h: b ≠ 0):
    emultiplicity p (a / b) = emultiplicity p a - emultiplicity p b := by
  admit -- defeat

-- omega instead of linarith
theorem imo2022_p5 {a b p : ℕ} (ha: 0 < a) (hb: 0 < b) (hp: Nat.Prime p): a^p = b ! + p ↔ (a, b, p) ∈ solutionSet := by
  refine ⟨fun h ↦  ?_, fun h ↦ ?_⟩
  . have gt_one : 1 < a := by
      apply lt_of_le_of_ne ha
      by_contra! hc
      replace hc : 1 = a := hc
      rewrite [← hc, one_pow] at h
      have := factorial_pos b
      have := hp.pos
      have : 1 < a := by omega
      exact this.ne.elim hc
    obtain hl | he | hg := lt_trichotomy a p
    . obtain hle | hgt := le_or_lt a b
      all_goals exfalso
      . have dvd_p : a ∣ a^p - b ! := dvd_sub' (dvd_pow_self a hp.ne_zero) $ dvd_factorial ha hle
        have : p = a^p - b ! := by omega
        rw [← this, dvd_prime hp] at dvd_p
        obtain rfl | ha' := dvd_p
        contradiction
        exact hl.ne ha'
      . have h' : b ! = a^p - p := by omega
        have lt_b_fac : a ! < a^p - p := by
          calc a !
            ≤ a^a - a := factorial_le_pow_sub gt_one -- fix
          _ < a^p - p := by sorry
        rw [← h'] at lt_b_fac
        rw [← factorial_lt hb] at hgt
        exact Nat.lt_asymm lt_b_fac hgt
    . rw [solutionSet]
      obtain rfl | ho := hp.eq_two_or_odd
      have : 2 ! = b ! := by simp_all
      rw [factorial_inj one_lt_two] at this
      rw [he, ← this]
      apply Set.mem_insert
      match p with
      | 2 => contradiction
      | 3 =>
        have : 4 ! = b ! := by simp_all; apply_rfl
        rw [he] at h ⊢
        rw [factorial_inj $ by norm_num] at this
        rw [← this]
        apply Set.mem_insert_of_mem
        rfl
      | 5 =>
        have : b ! = 3120 := by simp_all
        exfalso
        refine monotone_factorial.ne_of_lt_of_lt_nat 6 ?_ ?_ _ this <;> decide
      | p + 7 =>
        let p' := p + 7 -- change of variables
        have hp' : p + 7 = p' := by rfl
        rw [hp'] at ho ⊢
        have sub_one_even := ho.tsub_odd odd_one
        have mult_self: emultiplicity 2 2 = 1 := multiplicity.Finite.emultiplicity_self $ by decide
        have := Nat.two_pow_sub_pow sub_one_even.two_dvd ho.not_two_dvd_nat sub_one_even
        have : emultiplicity 2 (p' + 1) + emultiplicity 2 (p' - 1) + emultiplicity 2 (p' - 1) - 1 =
            emultiplicity 2 ((p' + 1) * (p' - 1) * ((p' - 1)/2)) := by
          calc
            _ = emultiplicity 2 ((p' + 1) * (p' - 1)) + emultiplicity 2 (p' - 1) - 1 := by rw [← emultiplicity_mul prime_two.prime]
            _ = emultiplicity 2 ((p' + 1) * (p' - 1)) + (emultiplicity 2 (p' - 1) - 1) := by admit
            _ = emultiplicity 2 ((p' + 1) * (p' - 1)) + (emultiplicity 2 (p' - 1) - emultiplicity 2 2) := by rw [mult_self]
            _ = emultiplicity 2 ((p' + 1) * (p' - 1)) + emultiplicity 2 ((p' - 1)/2) := by rw [emultiplicity_div two_ne_zero]
            _ = emultiplicity 2 ((p' + 1) * (p' - 1) * ((p' - 1)/2)) := by rw [← emultiplicity_mul prime_two.prime]

    -- Case 2: a > p
    . have le_pow : p ≤ p^p := by sorry -- Nat.le_self_pow hp.ne_zero
      have p_lt_a_pow := Nat.pow_lt_pow_left hg hp.ne_zero
      have p_lt_b : p ! < b ! := by
        calc p !
          ≤ p^p - p := factorial_le_pow_sub hp.two_le
        _ < a^p - p := Nat.sub_lt_sub_right le_pow p_lt_a_pow
        _ = b ! := by omega
      rw [Nat.factorial_lt] at p_lt_b
      have : p ∣ a^p := by
        rw [h]
        apply_rules [dvd_add, dvd_factorial hp.pos p_lt_b.le, dvd_rfl]
      have : p ∣ a := Nat.Prime.dvd_of_dvd_pow hp this
      sorry
      sorry

    sorry -- a ≠ 1
  rcases h with h0 | h1 | _
  simp_all
  apply_rfl
```

== Explanation

Tried to literally translate the paper proof into Lean code but ran into some implicit assumptions that were not stated clearly in the paper, such as $p! <= p^p - p$ or $a! < a^p - p$.