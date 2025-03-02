#import "notes.typ": *
#show: great-theorems-init

#show: notes.with(
  title: "Computer Architecture",
  subject: "Notes from the book",
  author: "Roland Yang"
)
#set heading(numbering: none)

#let blockquote(body) = {
  pad(left: 1em,
    block(
      stroke: (left: 1pt),
      inset: (left: 0.75em, y: 1em),
      body
    )
  )
}

#let enable_numbering = {
  set heading(numbering: "1.")
  show heading.where(level: 3): set heading(numbering: none)
}

= Introduction

These are a set of notes written for an independent study in Computer
Architecture led by Mr. Taylor Belcher at the Governor's School for Science and
Mathematics conducted in Spring 2025. We are using _Computer Architecture: A
Quantitative Approach (6th Edition)_.

From the book:

#blockquote[
  Brackets for each exercise (\<chapter.section>) indicate the text sections of
  primary relevance for completing the exercise. ... Exercises are rated, to
  give the reader a sense of the amount of time required to complete an exercise:

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

#enable_numbering
= Fundamentals of Quantitative Design and Analysis

== Introduction
- Rise of microprocessors resulted in
  - RISC architectures replacing others
  - Purpose-built computers being replaced
  - 50,000x performance boost traded language performance for productivity
- Dennard scaling: constant power density even as transistor density increased
  - Ended in 2004
- Cannot continue to rely on instruction-level parallelism (ILP)
- High performance uniprocessors #sym.arrow multiprocessors & data-level
  parallelism (DLP)
- _Amdahl's Law_ limits the performance benefit from parallelism as a function
  of how serial a task is
- Moore's law slowed from 1.5 years to 20 years for transistor doubling
- Only way to improve energy-performance-cost are specialized architectures

== Classes of Computers
- Five classes of computers:
  - IoT/Embedded: price over performance
  - Personal mobile device (PMD): cost, efficiency, limited by inadequate
    cooling, _real-time_ media requirements
  - Desktop: price-performance
  - Servers: availability, scalability, throughput
  - Clusers/Warehouse-Scale: price-performance, power, availability
- 2 kinds of application parallelism: _data-level parallelism_
  (DLP) and _task-level parallelism_ (TLP)
- Flynn (1996) named 4 types of parallelism
  - SISD
  - SIMD
  - MISD
  - MIMD

== Defining Computer Architecture

- "Old" view of computer architecture was simply ISA design
- "Genuine" architecture: ISA, _microarchitecture_, _hardware specifics (such
  as clock rate)_

== Trends in Technology
- To prepare for computer evolution, designers must follow changes in implementation tech
- Five are crucial to modern: _Integrated circuit logic_, _Semiconductor DRAM_,
  _Semiconductor Flash_, _Magnetic disk technology_, and _Network technolgy_
- Bandwith is improving faster than latency
- IC processes are characterized by _feature size_, the minimum size of
  transistor in the $x$ or $y$ direction.
- Shrinkage in the $y$ direction corresponds with reduction in operating voltage
  to maintain operation and reliability of transistors
- Transistor count improves quadratically with transistor performance
  - Density allowed for movement from 4-bit to eventually 64-bit computers,
    wide SIMD units, and innovations in speculative execution and caching
- Wire delay scales poorly compared to transistor performance

== Trends in Power and Energy in Integrated Circuits
- What is the maximum power a processor requires?
- What is sustained power consumption (TDP)?
  - _Thermal design power_: neither peak (1.5x) nor average—just enough for
    cooling
- Modern processors can manage heat by throttling and shutdown
- We compare processor efficiency via energy rather than power
  - Only use power consumption as a constraint of cooling, for example
- CMOS primary energy compution is switching transistors (dynamic energy)
  $ "Energy"_"dynamic" prop "Capacitive load" times "Voltage"^2 $
- Which is the energy of a transition of $0 arrow 1 arrow 0$; one transition
  is half
- Power per transistor is then transition energy times frequency
- Meaning slowing clock rate reduces power, not energy
- This relationship has meant that voltages have dropped from 5V to 1V in 20
  years
- Increase in transistor switching and higher frequency have outpaced decreases
  in load capacitance and voltage with newer processes, so power consumption has
  continued to grow
  - First microprocessors consumed $<1$W, while intel i6-6700K consumes 95W, we
    are at the limit of air cooling (and have been for a decade)
- Clock rate increase has slowed accordingly—how do we improve energy
  efficiency? We have four options:
  - Disabling clock of inactive cores or units
  - Dynamic voltage-frequency scaling (DVFS): reducing voltage when frequencies
    and voltages are unneeded
  - Designing for the typical case: think S sleep states where DRAM or disk
    powers off
  - Overclocking: running at higher clocks until tempuratures cause throttling,
    _Turbo mode_
- Though dynamic power is important, static power matters: leakage current flows
  even when off. $ "Power"_"static" prop "Current"_"static" times "Voltage" $
- Increasing transistor number increases power, and current leakage increases
  with transistor density
  - Resulted in turning off power supply entirely to submodules
- Use of faster, less-efficient processor to finish tasks and let system sleep
  faster (_race-to-halt_)
- Transistor density improvements have caused the creation of _dark silicon_,
  transistors that cannot be turned on because of thermal constraints
- Single precision FP add uses 30x the energy of 8-bit integer add
- Memory accesses (DRAM) use 20,000x more energy than 8-bit add
- Minimizing energy-per-task is up to domain specific processors

== Trends in Cost
- Costs go down over time as yields improve (the _learning curve_)
- Volume decreases time in learning curve and cost

== Dependability
- ICs were one of the most reliable components
  - As feature sizes shrink past 16 nm, transient and permanent faults grow more
    common
- Mean time to failure (MTTF) is a reliability measure
  - Reciprocal is _failures in time_ (FIT) in failures per billion hours
- Interruptions—_mean time to repair_ (MTTR)
- Mean time between failures (MBTF) is MTTF + MTTR
$ "Module availability" = "MTTF"/("MTTF" + "MTTR") $
- Cope with failure with redundancy (time or resources)

== Measuring, Reporting, and Summarizing Performance
