# Interpreting Benchmarks (version 0.1.7)

I'm using github style markdown here to present tables.

New benchmarks added for GF(2<sup>16</sup>) field.

On AMD A10:

| Test                |   Time    | Note                             |
| :---                |   ---:    | :---                             |
| gf16 mul/ref        | 144.84 us | Reference GF(2^16) mul           |
| gf16 mul/good       | 478.20 us | Much worse than reference        |
| gf16 inv/ref        | 14.833 us | Reference inv                    |
| gf16 inv/good       | 558.80 ns | Much better than reference       |
| gf16 div/ref        | 3.8748 ms | Reference GF(2^16) div           |
| gf16 div/good       | 599.82 us | Much better than reference       |


On Aarch64 (ODROID N2):

| Test                |   Time    | Note                             |
| :---                |   ---:    | :---                             |
| gf16 mul/ref        | 221.20 us | Reference GF(2^16) mul           |
| gf16 mul/good       | 1.1230 ms | Much worse than reference        |
| gf16 inv/ref        | 22.151 us | Reference inv                    |
| gf16 inv/good       | 1.5252 us | Much better than reference       |
| gf16 div/ref        | 5.9635 ms | Reference GF(2^16) div           |
| gf16 div/good       | 1.4071 ms | Much better than reference       |

## Assembly output (ARM)

Nothing unusual. I saw two bounds checks on mod reduce lookups (which
I eliminated) but none on MULL lookups. Probably because the compiler
determined that the lookup couldn't generate an invalid index. Nice.

## Analysis

With my C version of the same code, I do actually get faster results
with table lookups:

    Test 0: long (modular) multiply (u16 version)
    Test 0: -> 59.36284 M*/s
    Test 16: 8-bit x 4-bit straight *, l-r on nibbles, fast shift w/o init
    Test 16: -> 130.75171 M*/s

There are some differences, though:

* I use separate l-r tables
* I look up the l-r tables directly using indices constructed like `a3 | b1`

By using separate l-r tables, I would eliminate one shift on output
per rmull lookup, saving 4 instructions.

By using operands that are pre-shifted, with nibbles in the high four
bits, I would save a further 8 shifts.

Whether either of these will improve the situation enough to bring the
"good" mul implementation within striking distance of the reference
one isn't clear to me, though. I suspect it won't be enough.

Actually, looking at the generated assembly, it appears that the
compiler has already optimised the array indexing for me. 

### What's causing the discrepancy between Rust, C?

Since I can't see anything that's very inefficient in the table lookup
code, the other possibility is that the compiler has optimised the
*reference* code better than my C compiler has.

Right now, I don't have a way of looking at the assembly output
since the generic code isn't concretely instantiated in the library.

### Can other optimisations do any better than the reference code?

It should be possible. I chose this particular implementation of
GF(2<sup>16</sup>) to begin with because it has the smallest memory
footprint and because a similar technique can be applied to
GF(2<sup>32</sup>). Also, while a more obvious approach might have
been to use, say, full multiplication lookup tables or log/exp tables,
I'd already used those for GF(16) and GF(256), so I wanted to try
something different.

## Summary

More disappointing results, but at least I've added some more useful
building blocks.

As a high priority, I'll have to figure out how to emit assembly for
my generic field multiplication code.

I will leave both "good" multiplications where they are for the
moment, even though they are actually worse than the reference code.
When I have a more complete set of tables and I've refactored the code
to move them into `guff::tables`, I can set up more comprehensive
benchmarking of options and update `guff::good` according to the
results.
