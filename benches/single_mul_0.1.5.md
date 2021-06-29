# Interpreting Benchmarks (version 0.1.5)

I'm using github style markdown here to present tables.

Headline figures for AMD A10:

| Test                |   Time    | Note                             |
| :---                |   ---:    | :---                             |
| gf4 mul/ref         | 526.49 ns | Reference version of GF(16) mul  |
| gf4 mul/good        | 544.66 ns | Slightly worse than reference    |
| gf4 mul/mull        | 275.46 ns | Long multiply                    |
| gf4 mul/mull-reduce | 3.5463 us | Long mul + modular reduction     |
| gf4 inv/ref         | 153.82 ns | Reference version of inv         |
| gf4 inv/good        | 38.776 ns | Much better than reference       |
| gf8 mul/ref         | 112.05 us | Reference version of GF(256) mul |
| gf8 mul/good        | 137.36 us | Slightly worse than reference    |


On Aarch64 (ODROID N2):

| Test                |   Time    | Note                             |
| :---                |   ---:    | :---                             |
| gf4 mul/ref         | 1.2899 us | Reference version of GF(16) mul  |
| gf4 mul/good        | 1.2116 us | Slightly better than reference   |
| gf4 mul/mull        | 711.73 ns |                                  |
| gf4 mul/mull-reduce | 6.2504 us |                                  |
| gf4 inv/ref         | 245.31 ns | Reference version of inv         |
| gf4 inv/good        | 64.380 ns | Much better than reference       |
| gf8 mul/ref         | 184.77 us | Reference version of GF(256) mul |
| gf8 mul/good        | 294.09 us | Much worse than reference        |

## Assembly output (ARM)

Using `cargo asm` to examine the GF(2<sup>8</sup>) multiply code, I
see nothing unusual:

```Assembly
	<guff::good::F8_0x11b as guff::GaloisField>::mul:
	ldr     x8, [x0]               // addr of log table 
	and     x9, x1, #0xff          // operand a
	and     x10, x2, #0xff         // operand b
	ldrsh   x9, [x8, x9, lsl, #1]  // log[a]
	ldrsh   x8, [x8, x10, lsl, #1] // log[b]
	ldr     x10, [x0, #48]         // addr of exp_entry 
	add     x8, x8, x9             // log a + log b 
	ldrb    w0, [x10, x8]          // exp[sum]
	ret
```

The important thing to note here is that there are 5 memory accesses.
Two are for getting addresses of the `exp` and `log` tables, and then
three for looking up values in those tables.

There's also a bit of overhead applying a mask of `0xff` to the
operands. I'm not sure why that should be, since the inputs are of
type `u8`. 

The assembly output on `x86_64` is similar.

On both platforms, it seems that the costs of memory accesses is
(broadly speaking) higher than the cost of calculating the product
using bit manipulation and loop counting. Also, the extra masking that
happens in the ARM code is probably contributing somewhat to worse
performance.

## Possible improvement: static tables

We can remove two memory lookups by storing static tables within the
library. I will implement this for the fields in `guff::good`, but I
want to avoid cluttering up (bloating) the library with lookup tables
that probably won't be used. I can use Rust's `feature` mechanism to
make static tables optional.

I'll need to do a bit of refactoring to make that work. Another,
related bit of refactoring will be moving the particular
implementation of fields out of `guff::good` and into a new module
`guff::tables`. 

As an interim measure, so as to avoid needing to refactor anything too
soon, I will start adding static tables into the `guff::good`
module. I won't include the code to generate those tables at first,
though. In the longer term, though, I will:

* write Rust code for generating the tables (port from C code I have
  already written)

* move static tables into a new module and allow them to be selected
  as features

The `guff::good` module will continue to act as a kind of switchboard,
plugging into whatever static tables are available, or falling back to
dynamically-generated tables if not.

## It's not all bad

I was a little bit surprised that the "accelerated" multiply routines
were, in most cases, worse than the reference implementation. However,
I'm not too concerned about it right now.

For one thing, we probably shouldn't be too focused on how fast we can
do single multiply operations. It's fine for benchmarking, but in
reality, most users will want to work with vectors and matrices. In
that case, looking up table addresses (from some memory location) will
only need to be done once, at the start of each vector or matrix
operation. When this is amortised over length of the vector or size of
the matrix, it will mostly become negligible.

The second reason for not being too worried is that for smaller field
sizes, it's easy for any potential gains with table lookups to get
wiped out by the practicalities of memory accesses. When I move on to
implementing table lookups for GF(2<sup>16</sup>) and
GF(2<sup>32</sup>), I expect that the "good" implementations will
start to pull ahead.

Thirdly, and finally, my ultimate destination with this library is to
include SIMD code for working on vectors and matrices. I've already
implemented some of this in assembly, and the results are quite
impressive. The ARM NEON code running on a moderately powerful armv7
board (ODROID XU4) already outperforms my best software-based (single)
multiply routine running on my AMD A10 desktop. Though obviously, I
shouldn't read too much into that, since, as I've said above,
benchmarking single multiply throughput isn't really that meaningful.

## Updates

I will probably add more documentation like this from time to time. I
won't necessarily do it every release, though.

Where I do include benchmarks, I will use the official release code so
that others can reproduce them.
