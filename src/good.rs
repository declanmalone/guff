//! # A set of "good" optimised routines in native Rust
//!
//! This module provides a set of alternative implementations of
//! [GaloisField] that are faster than the reference implementation.
//! This is provided for those who are interested in having a
//! reasonable set of fairly fast routines, but who don't particularly
//! care about how they're implemented.
//!
//! This module is intended to provide a reasonable balance between:
//!
//! * speed;
//! * memory footprint;
//! * range of field sizes; and
//! * supported polynomials
//!
//! The actual implementation used here may change over time. In fact,
//! some of the fields may even default to using the reference
//! implementation for some calculations. It can be assumed, however,
//! that at least the `mul` routines implemented by objects here will
//! be faster than the corresponding reference implementation, since
//! this is usually the method that will be called most often in an
//! application.
//!
//! # SemVer implications
//!
//! Once I have added a constructor here, it will continue to be
//! available across all minor versions of the library. In other
//! words, in terms of semantic versioning, removing a constructor
//! here counts as breaking the application interface.
//!
//! However, no guarantees are made as to how the optimised version
//! works internally. In particular, no guarantee of speed is made.
//! Nor should the contents of the `struct` implementing the method be
//! considered part of the application interface. That can change at
//! any time.
//!
//! # Benchmarking
//!
//! I've implemented some basic benchmarks using [Criterion]. After
//! downloading the source for this project, simply run:
//!
//!```ascii
//!     cargo bench
//!```
//!
//! Tests marked `ref` are run using the default, reference
//! implementations, while those marked with `good` are for fields
//! constructed using the functions below.
//!

use crate::{ GaloisField };
use num::{One,Zero};
//use num_traits;
//use num_traits::ToPrimitive;

// For the first version, I'm going to use a variety of different
// techniques for each field size. I intend to migrate the code for
// calculating tables into a separate module, so if I use different
// techniques for each field, that means I'll have more code written
// and tested, ready for the migration.
//
// GF(2<sup>4</sup>):
//
// * Full mul tables for `mul` (256 bytes)
// * full inverse table for `inv`
// * rest supplied by default
//


// We appear to have to jump through some hoops to make code generic
// wrt GaloisField. However, it works!
struct FullMulLUT<G> where G : GaloisField {
    table  : Vec<G::E>,
}
impl<G> FullMulLUT<G>
where G : GaloisField, G::E : Into<usize>
{
    fn new(f : &G) -> FullMulLUT<G> {
	let zero = G::E::zero();
	let one  = G::E::one();
	let max : usize = (1 << (f.order() as usize)) - 1;
	let mut v = Vec::<G::E>::with_capacity((max + 1) * (max + 1));

	// let a correspond to rows; fill matrix rowwise
	let mut a = zero;
	let mut b;

	// eprintln!("order is {}", f.order());
	// eprintln!("max is {}", max);

	for _row in 0..=max {
	    b = zero;
	    for _col in 0..=max {
		let prod = f.mul(a,b);
		// eprintln!("a: {}, b: {}, prod: {}", a, b, prod);
		v.push(prod);
		b = b + one;
	    }
	    a = a + one;
	}
	FullMulLUT::<G> { table : v }
    }
    #[inline(always)]
    fn mul(&self, a : G::E, b : G::E) -> G::E {
	let index : usize = (a.into() << G::ORDER) + b.into();
	// Can't do unsafe access without ensuring a,b < 16: 
	// unsafe {
	//    *self.table.get_unchecked(index)
	// }
	self.table[index]
    }
}

// Note that I didn't have to make the above generic on a particular
// GaloisField implementor, and in fact it probably did make the job
// harder than it should have been. Here's an alternative way. It
// doesn't cut our dependence on GaloisField, but at least there's
// less faffing around getting max (and other counting numbers, if we
// had needed them)

// Also, note that we're just using a function instead of a full
// struct/impl combo.

fn fill_inverse<T>(f : & T,
		   v : &mut Vec<T::E>, max : usize)
    where T : GaloisField
{
    // eprintln!("max is {}", max);
    let mut elem = T::E::zero();
    v.push(elem);
    for _count in 1..=max {
	elem = elem + T::E::one();
	v.push(f.inv(elem));
    }
}

// "good" F4 with fixed poly 0x13 using above mul table
// Not meant to be used directly. Use [new_gf4_0x13] constructor
// instead.
#[doc(hidden)]
pub struct F4_0x13 {
    // Note how we're treating F4 solely as a type
    // (it's never allocated or stuffed in our struct)
    mul_lut : FullMulLUT::<crate::F4>,

    // Alternative method of construction used for inv_lut:
    inv_lut : Vec<u8>,

    // don't need to store poly (put in impl method calls instead)
}

impl GaloisField for F4_0x13 {
    type E = u8;
    type EE = u8;
    type SEE = i8;

    // we have to redeclare types for constants
    const ORDER      : u16 = 4;
    const POLY_BIT   : u8  = 0x10;
    const FIELD_MASK : u8  = 0x0f;
    const HIGH_BIT   : u8  = 0x08;

    // the two required methods (everything else is default)
    fn poly(&self)      -> u8  { 3 }
    fn full_poly(&self) -> u8  { 19 }

    // pass mul call on to table
    fn mul(&self, a : Self::E, b : Self::E) -> Self::E {
	self.mul_lut.mul(a,b)
    }

    fn inv(&self, a : Self::E) -> Self::E
    {
	// can use 'a as usize' since its type is known to be u8
	self.inv_lut[a as usize]
    }
}

/// Optimised maths for GF(2<sup>4</sup>) with the (primitive) polynomial 0x13
pub fn new_gf4_0x13() -> F4_0x13 {
    // reference field object
    let f = crate::new_gf4(19,3);
    // generate inverse table
    let mut inv = Vec::<u8>::with_capacity(16);

    fill_inverse(&f, &mut inv, 15);
    
    F4_0x13 {
	mul_lut : FullMulLUT::<crate::F4>::new(&f),
	inv_lut : inv,
    }
}




// GF(2<sup>8</sup>):
//
// * extended log/exp tables for `mul`,`div`,`inv`, `pow`
// * optimise multiplying vector by constant
// * rest supplied by default

// This presents a difficulty because we don't store signed versions
// of types in GaloisField. It's not beyond the bounds of possibility
// that we would use log/exp tables for field sizes up to u16. As a
// result, it's probably worth making the table creation function
// generic instead of repeating it for u4, u8 and u16.

struct BigLogExpTables<G> where G : GaloisField {
    log  : Vec<G::SEE>,
    exp  : Vec<G::E>,
    exp_entry : *const G::E,
}
impl<G> BigLogExpTables<G>
where G : GaloisField,
      G::E : Into<usize>,
      G::E : Into<G::SEE>,
      G::SEE : Into<isize>,
//      G::SEE : From<isize>,
//      G::SEE : From<usize>,
      G::E : std::fmt::Debug
{
    fn new(f : &G, g : G::E ) -> BigLogExpTables<G> {

	// eg, for GF256, log_size = 256, exp_size = 1024
	let log_size = 1 << (G::ORDER as usize);
	let see_log_size : G::SEE = G::SEE::one() << (G::ORDER as usize);
	let exp_size = log_size * 4;
	let exp_entry;

	// add to exp sequentially, to log randomly
	let zero = G::SEE::zero();
	let mut exp : Vec<G::E>   = Vec::<_>::with_capacity(exp_size);
	let mut log : Vec<G::SEE> = vec![zero ; log_size];

	// for gf(256), add 512 zeros
	for _ in 0..log_size * 2 {
	    exp.push(G::E::zero());
	}
	unsafe {
	    // offset 512 from start of exp, which is within bounds
	    exp_entry = exp.as_ptr().offset(log_size as isize * 2);
	}

	// first bunch of entries
	let mut i : usize = 0;
	// exp 0 = 1
	exp.push(G::E::one());
	// log 0 = -256
//	log[i] = (- (log_size as isize)).into();
	log[i] = G::SEE::zero() - see_log_size;

	// exp 1 = generator
	exp.push(g);
	// log g = 1
	let usize_g : usize = g.into();
	log[i + usize_g] = G::SEE::one();

	// usize loop counter and G::SEE one
	i += 2;
	let mut ei = G::SEE::one() + G::SEE::one();
	let mut p = g;			// running product
	loop {
	    if p == G::E::one() {
		panic!("{} is not a generator for this field", g)
	    }

	    p = f.mul(p,g);
	    exp.push(p);
	    let usize_p : usize = p.into();
	    log[usize_p] = ei;
	    i += 1; ei = ei + G::SEE::one();
	    if i == log_size { break }
	}

	// We have inserted all log table entries.
	// 
	// exp entries currently stand at:
	//
	// 512 zero values
	// g ** 0 = 1
	// g ** 1 = g
	// ...
	// g ** 255 = 1
	//
	// We now have to add values starting at g again, finally
	// finishing with another g (last two entries are optional:
	// since max sum of logs is 255 + 255, they will never be
	// accessed from mul)
	assert_eq!(p, G::E::one());
	for _ in 0..log_size-1 { // 
	    p = f.mul(p,g);
	    exp.push(p);
	}
	assert_eq!(p, G::E::one());
	// last element must be zero for inv(0) = 0 to work
	exp.push( G::E::zero() );
	assert_eq!(exp_size, exp.len());
	
	BigLogExpTables::<G> { log, exp, exp_entry }
    }
    #[inline(always)]
    fn mul(&self, a : G::E, b : G::E) -> G::E
    {
	let usize_a : usize = a.into();
	let usize_b : usize = b.into();
	let log_a : isize;
	let log_b : isize;
	unsafe {
	    // safe because log table has entry for each field element
	    log_a = (*self.log.get_unchecked(usize_a)).into();
	    log_b = (*self.log.get_unchecked(usize_b)).into();
	    // safe because log_a + log_b within exp table bounds:
	    // -512 ... 510
	    *(self.exp_entry.offset(log_a + log_b))
	}
	// replace with unsafe after testing ...
	//	self.exp[(512 + log_a + log_b) as usize]
    }
    // can also implement inv, div, pow with these tables!
    fn inv(&self, a : G::E) -> G::E {
	// for GF(256), we would return exp[255 - log a]
	let log_top = (1 << (G::ORDER as usize)) - 1;
	let usize_a : usize = a.into();
	let log_a : isize;
	unsafe {
	    log_a = (*self.log.get_unchecked(usize_a)).into();
	    *(self.exp_entry.offset(log_top - log_a))
	}
    }
}

// I will implement two fields here:
//
// 0x11b : non-primitive. It's the one I often use. It's also the poly
// used in AES. Has generator `3`
//
// 0x11d : primitive poly, so could be more useful in other
// applications. Since it's primitive, the generator is `2`

// I'll follow the same process as with the GF(16) poly above. Note
// that no memory is allocated for any of these fields unless the user
// calls the constructor. There is the overhead for making a concrete
// u8 version of the BigLogExpTables code, though.

#[doc(hidden)]
pub struct F8_0x11b {
    tables : BigLogExpTables::<crate::F8>,
}

impl GaloisField for F8_0x11b {
    type E = u8;
    type EE = u16;
    type SEE = i16;

    // we have to redeclare types for constants
    const ORDER      : u16 = 8;
    const POLY_BIT   : u16 = 0x100;
    const FIELD_MASK : u8  = 0xff;
    const HIGH_BIT   : u8  = 0x80;

    // the two required methods (everything else is default)
    fn poly(&self)      -> u8   { 0x1b }
    fn full_poly(&self) -> u16  { 0x11b }

    // pass mul call on to table
    // #[inline(always)] // doesn't play nicely with benchmark code
    fn mul(&self, a : Self::E, b : Self::E) -> Self::E {
	self.tables.mul(a,b)
    }

    fn inv(&self, a : Self::E) -> Self::E
    {
	self.tables.inv(a)
    }
}

/// Optimised maths for GF(2<sup>8</sup>) with the (non-primitive) polynomial 0x11b
pub fn new_gf8_0x11b() -> F8_0x11b {
    // reference field object
    let f = crate::new_gf8(0x11b,0x1b);
    
    let this = F8_0x11b {	// field has generator 3
	tables : BigLogExpTables::<crate::F8>::new(&f, 3),
    };

    // Ensure that log[0] == -256
    assert_eq!(this.tables.log[0], -256);

    this
}



//
// GF(2<sup>16</sup>):
//
// * l-r with 8-bit modular shift, breaking operands into four nibbles
//   and two bytes for `mul`
// * rest supplied by default
//
// GF(2<sup>32</sup>):
//
// * as per 16-bit, but breaking both operands into four 8-bit values
//   for `mul`
// * rest supplied by default
//
// 

// GF(2<sup>4</sup>) field implementations
//
// There are only three polynomials of order 4 that are irreducible
// for this field:
//
// * 0x13 (19) (primitive)
// * 0x19 (25) (primitive)
// * 0x1f (31)
//
// There's probably no reason to use 0x1f. 

#[cfg(test)]
mod tests {

    use super::*;
    use crate::{new_gf4, new_gf8};

    #[test]
    fn test_f4_0x13_mul_conformance() {
	let f4      = new_gf4(19,3);
	let f4_0x13 = new_gf4_0x13();
	let mut fails = 0;
	for i in 0..16 {
	    for j in 0..16 {
		if f4.mul(i,j) != f4_0x13.mul(i,j) {
		    fails += 1;
		}
	    }
	}
	assert_eq!(fails, 0);
    }

    #[test]
    fn test_f4_0x13_inv_conformance() {
	let f4      = new_gf4(19,3);
	let f4_0x13 = new_gf4_0x13();
	let mut fails = 0;
	for i in 0..16 {
	    if f4.inv(i) != f4_0x13.inv(i) {
		fails += 1;
	    }
	}
	assert_eq!(fails, 0);
    }

    #[test]
    fn test_f8_0x11b_mul_conformance() {
	let f8       = new_gf8(0x11b,0x1b);
	let f8_0x11b = new_gf8_0x11b();
	let mut fails = 0;
	for i in 0..=255 {
	    for j in 0..=255 {
		if f8.mul(i,j) != f8_0x11b.mul(i,j) {
		    fails += 1;
		}
	    }
	}
	assert_eq!(fails, 0);
    }

    #[test]
    fn test_f8_0x11b_inv_conformance() {
	let f8       = new_gf8(0x11b,0x1b);
	let f8_0x11b = new_gf8_0x11b();
	let mut fails = 0;
	for i in 0..=255 {
	    if f8.inv(i) != f8_0x11b.inv(i) {
		eprintln!("Failed inv({})", i);
		assert_eq!(f8.inv(i), f8_0x11b.inv(i), "(ref vs good");
		fails += 1;
	    }
	}
	assert_eq!(fails, 0);
    }
}
