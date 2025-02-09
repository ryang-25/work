#import "notes.typ": *
#show: great-theorems-init

#let blockquote(body) = {
  pad(left: 1em,
    block(
      stroke: (left: 1pt),
      inset: (left: 0.75em, y: 1em),
      body
    )
  )
}

#show: notes.with(
  title: "Computer Architecture",
  subject: "Notes from the book",
  author: "Roland Yang"
)
#set heading(numbering: none)

= Introduction

These are a set of notes written for an independent study in Computer Architecture led by Mr. Taylor Belcher at the Governor's School for Science and Mathematics conducted in Spring 2025. We are using _Computer Architecture: A Quantitative Approach (6th Edition)_.

From the book:

#blockquote[
  Brackets for each exercise (\<chapter.section>) indicate the text sections of primary relevance for completing the exercise. ... Exercises are rated, to give the reader a sense of the amount of time required to complete an exercise:

  [10] Less than 5 min (to read and understand) \
  [15] 5-15 min for a full answer \
  [20] 15-20 min for a full answer \
  [25] 1 h for a full written answer \
  [30] Short programming project: less than 1 full day of programming \
  [40] Significant programming project: 2 weeks of elapsed time \
  [Discussion] Topic for discussion with others
]

Our initial plan of progression is as follows:

#align(center)[Appendix A $arrow.r$ Chapter 1 $arrow.r$ Appendix B $arrow.r$ Chapter 2 $arrow.r$
Chapter 3 $arrow.r$ Chapters 4-7]

with additional appendices, as appropriate. Each corresponding portion of Appendix M will be read after completing each chapter.


= Appendix A

== 1.

- Desktop computing cares about integer and floating point performance, not code size
- Servers focus on integers and character strings rather than FP performance
- Personal and embedded care about code size because of memory costs and may forgo FP isns or use a compressed ISA
- Intel chips use a RISC ISA internally for performance while still presenting x86, but there are disadvantages to doing so, discussed later

== 2.

- 3 types of ISA:
  - Stack
  - Accumulator
  - General Purpose Registers (GPR)
- Implicit v explicit operands
- Stack: implicitly on top of stack
- Accumulator: one operand implicitly accumulator
- GPR: only explicit (registers or memory)
- Explicit operands may either be from memory or loaded into storage
- 2 classes of register computers:
  - Register-memory: can access memory w/ any isn
  - Load-store: can only access memory w/ loads and stores
- Secret third class (memory-memory)
  - All operands in memory
- Some ISA have more registers than accumulator but restricted
  - _extended accumulator_ or _special-purpose register_
- Early computers used stack or accumulator-style ISA
  - New architectures (>1980) used load-store
  - Why? Registers faster than memory
    - More efficient for compilers
    - `(A*B) + (B*C) - (A*D)` can be evaluated in any order in registers
    - On a stack computer can only be evaluated in one order
  - Registers can hold variables
    - When registers hold variables:
      - Memory traffic is reduced
      - Program speeds up
      - Code density is improved
  - Modern compilers have led to increase in register counts
  - Two instruction set characteristics divide GPR ISAs
    - Whether an ALU instruction has 2 or 3 operands
      - In 2 operands, one operand is both source and result
    - How many operands may be memory addresses in ALU isns (ranges from 0-3)

== 3.

- Little endian vs big endian
- Aligned accesses are faster than non-aligned

== 4.

== 5.

- Most widely executed isns are the simple ops

== 6.

- 4 types of control flow change
  - Conditional branches
  - Jumps
  - Calls
  - Returns
- Destination addresses for control flow isns are explicit
- Most common way are PC-relative
  - Allows PIC
- If not PC-relative need dynamic addressing
  - Register indirect jumps or other addressing mode
- Register indirect jumps useful for case/switch statements, virtual calls, function pointers, and dynamically shared libraries
- Many ISAs have compare and branch because of ubiquity

== 7-12

...

== Selected Exercises: 1, 7, 10


=== 1. Effective CPI

`gcc`: $.17*5 + .23*3 + .2*4 + .04*3 + .36*1 = 2.82$ CPI \
`astar`: $.27*5 + .06*3 + .18*4 + .02*3 + .46*1 = 2.77$ CPI

Loads Stores Branches Jumps ALU

= Fundamentals of Quantitative Design and Analysis

== 1.

- Rise of microprocessors resulted in RISC architectures replacing other machines
- Purpose-built computers replaced by microprocessors
- 50000x performance boost allowed programmers to trade performance for productivity
- Dennard scaling: constant power density even as transistor density increased
  - Ended in 2004
- High performance uniprocessor projects #sym.arrow multiprocessors and data-level parallelism
- Moore's law slowed from 1.5 years to 20 years for performance doubling

== 2.

- 2 kinds of application parallelism: _data-level parallelism_
  (DLP) and _task-level parallelism_ (TLP)
- Flynn (1996) named 4 types of parallelism
  - SISD
  - SIMD
  - MISD
  - MIMD

== 3.


