# Grand Unified Finite Field library*

Implements GF(2<sup>x</sup>) for various "natural" sizes
such as 2<sup>8</sup>, 2<sup>16</sup> and 2<sup>32</sup>.

My goals for this crate are to:

1. help me learn to write good modules in Rust;

2. help interested users learn about finite fields (ie, Galois
fields);

3. provide a generic baseline implementation of basic maths
(add, multiply, divide, etc.) over finite fields;

4. explore various optimisations/adaptations (including
table-based lookups and architecture-specific SIMD code) that can
selectively override some/all of the default implementations
while remaining compatible with other implementations.

Also to:

5. provide some useful utility functions that go beyond just
`add`, `mul`, `div`, etc. (eg, determining whether a field
polynomial is primitive, or generating lookup tables for different
kinds of optimisations);

6. eat my own dog food, so to speak, by implementing various
   applications that use the library;

7. benchmark particular implementations of interest.

# Basic Use: doing maths in a particular field

Steps to using this library:

* decide what "class" of field you want to use (GF(2<sup>8</sup>),
GF(2<sup>16</sup>), etc.);

* decide if you want to use one of the optimised adaptations or
are happy with the default generic code;

* create a new field object (we can call `f`) of that class with
your chosen field polynomial (aka "irreducible polynomial") by
calling the appropriate constructor;

* use that object to do maths in that field: eg, `result =
f.mul(a,b)`

# Future Work

## Vectors

- [ ] Decide on function suite and give better names
- [ ] Write test cases
- [ ] Write benchmarks

## SIMD

- [x] Write and test x86 Rust code (using `core::arch`)
- [ ] Import that code into the project
- [ ] Use CPU/feature detection to conditionally compile it
- [x] Write and test Arm NEON code in C
- [ ] Port that code to Rust
- [ ] Implement fallback code for unsupported platforms

## Refactoring

- [ ] Move table generation into `guff::tables`
- [ ] Encapsulate inverse table code in a similar style
- [ ] New module `guff::impls` for different ways to implement field maths
- [ ] Use preamble (?) for list of sensible exports
- [ ] Make some large static tables available as selectable features

## Extra functionality

- [ ] Finalise set of table generation routines
- [ ] Test whether a value is a generator for a field
- [ ] Test whether a polynomial is irreducible
- [ ] Test whether a polynomial is primitive
- [ ] Reference lists of polynomials

## Pure Rust "good" implementations

- [ ] Benchmark full set of optimised implementations
- [ ] Update `guff::good` with results

## Documentation

- [ ] Write short intro to finite fields, with bibliography
- [ ] Gather list of similar projects and write summaries/mini reviews
- [ ] Write about experience of writing this project in Rust
- [ ] Write short article about Assembly/SIMD in Rust

# Crate Name

\* The crate name is deliberately hyperbolic:

> Noun *guff* - unacceptable behavior (especially ludicrously false statements)

# Copyright and Licence

This work is Copyright (c) Declan Malone, 2021.

You may freely copy and modify this work under the terms of:

* The GNU General Public License version 2 or later

If you wish to embed this work as part of another work, you may do so
under the terms of:

* The GNU Lesser General Public License version 2 or later

Disclaimer: this software comes with no warranty, express or implied.
