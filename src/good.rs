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
//! implementation for some calculations. 
//!
//!
//!
//!
//!
//!
//!
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
    bits   : usize,
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
	FullMulLUT::<G> { bits : f.order() as usize, table : v }
    }
    fn mul(&self, a : G::E, b : G::E) -> G::E {
	let index : usize = (a.into() << self.bits) + b.into();
	self.table[index]
    }
}

// Note that I didn't have to make the above generic on a particular
// GaloisField implementor, and in fact it probably did make the job
// harder than it should have been. Here's an alternative way. It
// doesn't cut our dependence on GaloisField, but at least there's
// less faffing around getting max (and other counting numbers, if we
// had needed them)

fn fill_inverse<T>(f : & T,
		   v : &mut Vec<T::E>, max : usize)
    where T : GaloisField
{
    eprintln!("max is {}", max);
    let mut elem = T::E::zero();
    v.push(elem);
    for _count in 1..=max {
	eprintln!("_count: {} ", _count);
	elem = elem + T::E::one();
	eprintln!("_count: {}, elem: {}", _count, elem);
	v.push(f.inv(elem));
	eprintln!("_count: {} ", _count);
    }
	panic!();
}

// "good" F4 with fixed poly 0x13 using above mul table
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
    use crate::new_gf4;

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

}
