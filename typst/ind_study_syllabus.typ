// formatting tweaks

#set page(paper: "us-letter", margin: 1in)
#set text(11pt)
#show heading.where(level: 1): text.with(size: 12pt)
#show heading: text.with(font: "Inter", size: 10pt)
#show heading: set block(above: 1.25em, below: 1em)

#set par(spacing: 1em)
#set list(indent: .75em)
#set enum(indent: .75em)

// header
#image("assets/logo.png", alt: "GSSM Logo", width: 30%)

= Course Syllabus

Seminar in Advanced Computer Science (Independent Study) \
Spring 2025

#grid(
    columns: (1.25in, auto),
    row-gutter: 0.75em,
    [Instructors:], [Dr. Al DeGennaro and Mr. Reid Mewborne],
    [Email:], [#link("mailto:degennaro@governors.school")],
    [Phone:], [#link("tel:8433833891,1119")[(843) 383-3901 x1119]],
    [Office:], [C119],
    [Zoom ID:], [#link("https://gssm.zoom.us/j/4335620678")[433-562-0678]]    
)

_Prerequisites:_ Completion of *CSC230* (Data Structures and Algorithms)

== Course Overview

This is an independent study course meeting weekly over videoconference covering
topics in low-level computer science, including computer organization, operating
system design, and hardware interfacing. Students are expected to have strong
theoretical and practical computer science backgrounds, including:

- Experience with low-level application programming in a language such as C/C++
- Experience debugging their code with portable debuggers such as gdb or lldb
- Experience writing assembly language for 64-bit architectures
- Awareness of data structures, algorithm design, and computational complexity

Among other tasks, students will fully design and implement a modern operating
system (OS) from scratch in #link("https://www.rust-lang.org")[Rust], a new
systems programming language, making use of its strong memory safety guarantees
to write safe, performant code. Students will be expected to follow recent
advances in extant operating systems during implementation and provide strong
empirical rationale—through benchmarking—for their design decisions. Their OS
will target a commodity single-board computer, the Raspberry Pi 4, running on
the AArch64 architecture.

During meetings, students should expect to discuss their readings and progress
on current projects with instructors, as appropriate. Through discussion and
project-based work, students will leave with a greater understanding in computer
organization and other advanced topics in computer science.

== Course Objectives

At the end of the course, students will be able to:

+ Interact and navigate the git version control system and command line utility
  to version and maintain copies of code.
+ Effectively write and review codebases and changes, accepting and offering
  valuable feedback in corporate settings.
+ Be proficient in low-level programming in the Rust programming language and
  its ownership and borrowing primitives as a tool for writing memory-safe code.
+ Utilize linker script to craft binaries suitable for booting on single-board
  computers.
+ Create a fully functional, minimal, extensible operating system.
+ Create and present a final presentation summarizing their work.

#grid(
  columns: (1.3cm, auto),
  row-gutter: 0.75em,
  [*Texts:*],
  [_The Rustonomicon_., found at
  #link("https://doc.rust-lang.org/nomicon/")[doc.rust-lang.org/nomicon/].],
  [], [_Arm Architecture Reference Manual for A-profile architecture_,
  #link("https://developer.arm.com/documentation/ddi0487/latest")[developer.arm.com/documentation/ddi0487/latest].],
  [], [Instructors will provide other texts as necessary.]
)


== Course of Study

May be subject to change as the discretion of the instructors.

#table(
  columns: (auto, 1fr),
  [Week 1], [ARM Architectural Timers],
  [Week 2], [JTAG Hardware Debugging],
  [Week 3], [Exception Levels (EL)],
  [Week 4], [Memory Management Unit (MMU) / Virtual Memory],
  [Week 5], [CPU Exceptions],
  [Week 6], [Unit and Integration Testing],
  [Week 7], [Asynchronous Exceptions],
  [Week 8], [Userland Virtual Memory],
  [Week 9], [Translation Tables],
  [Week 10], [Higher-Half Kernel],
  [Week 11], [Symbols],
  [Week 12], [Backtracing],
  [Week 13], [Dynamic Memory Allocation],
  [Week 14], [Timer Callbacks],
  [Week 15], [TBD],
  [Week 16], [Final Project Preparation]
)
