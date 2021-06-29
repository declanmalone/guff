//! # Grand Unified Finite Field library*
//!
//! Implements GF(2<sup>x</sup>) for various "natural" sizes
//! such as 2<sup>8</sup>, 2<sup>16</sup> and 2<sup>32</sup>.
//! 
//! # Basic Use: doing maths on elements of a particular field
//! 
//! The general outline for using this library is:
//!
//! * decide what "class" of field you want to use (GF(2<sup>8</sup>),
//! GF(2<sup>16</sup>), …);
//! 
//! * decide if you want to use one of the optimised adaptations or
//! are happy with the default, generic code;
//!
//! * create a new field object (we can call `f`) of that class with
//! your chosen field polynomial (aka "irreducible polynomial") by
//! calling the appropriate constructor;
//!
//! * use that object to do maths in that field: eg, `result =
//! f.mul(a,b)`
//! 
//! Example of using the default generic code:
//! 
//! ```rust
//! use guff::{GaloisField, F4, new_gf4};
//! 
//! // Create a GF(2<sup>4</sup>) field from a struct
//! // * `19` is the field's irreducible polynomial
//! // * `3` is the same value with bit 0x10 removed
//! let f = guff::F4 { full : 19, compact : 3 };
//! 
//! assert_eq!(f.pow(5,3), f.mul(5,f.mul(5,5)) );
//! 
//! // The same, but using a constructor as syntactic sugar
//! let f2 = guff::new_gf4(19, 3);
//!
//! assert_eq!(f2.pow(5,3), f2.mul(5,f2.mul(5,5)) );
//! assert_eq!(f.pow(5,3), f2.pow(5,3));
//! 
//! ```
//! 
//! The `guff::good` module provides a set of reasonably good
//! implementations that in most cases *should* be better than the
//! generic implementation. Example of use:
//! 
//! ```rust
//! use guff::GaloisField;
//! use guff::good::{F4_0x13, new_gf4_0x13};
//! 
//! let f = guff::good::new_gf4_0x13();
//!
//! assert_eq!(f.pow(5,3), f.mul(5,f.mul(5,5)) );
//! 
//! ```
//!

//! 
//! # Vector Operations
//! 
//! Many applications involving Galois Fields involve working with
//! vectors of field elements instead of single elements. Various
//! primitive methods are implemented:
//! 
//! * sum of vector elements
//! * sum of two equal-length vectors
//! * cross product (pairwise multiplication of elements)
//! * dot product (sum of cross-product values)
//! * scaling a vector by a constant
//! * fused multiply add (scale and sum by a pair of constants across
//!   vector)
//! 
//! These are all implemented using slices of the appropriate
//! [ElementStore] type. Where necessary, if a vector (slice) type is
//! to be returned, it must be handled by the user by passing in a
//! mutable reference.  Also note that some operations (such as
//! scaling a vector) are done in-place, so if the previous values
//! need to be saved, they should be copied first.
//! 
//! See the [GaloisField] documentation for a full list of method
//! signatures (prefixed with `vec_`).
//! 
//! # Crate Name
//! 
//! \* The crate name is deliberately hyperbolic:
//!
//! > Noun *guff* - unacceptable behavior (especially ludicrously
//!   false statements)

//use num_traits;
use num::{PrimInt,One,Zero};

pub mod good;


// I think that if I want to keep a flat directory structure, while
// still supporting an arbitrarily deep module tree, I would have to
// create `tables.rs` and have a `pub mod mull` line within it.

// put mull into tables
pub mod mull;      // I thought I could make this private? No?
pub mod tables {
    pub use crate::mull;
}

// pub use mull as tables::mull;

/// A typing trait meant to map to a primitive unsigned integer type
/// such as u8, u16 or u32.
pub trait ElementStore : 'static + Copy
    + num_traits::int::PrimInt
    + std::fmt::Display
    + num::FromPrimitive + num::ToPrimitive
    // + num::FromPrimitive + num::ToPrimitive
    // + num::Integer
    // + num::cast::AsPrimitive<Self::E>
    // + num::Integer + num::One
{}
impl<T> ElementStore for T
where T : 'static + Copy
    + num_traits::int::PrimInt
    + std::fmt::Display
    + num::FromPrimitive + num::ToPrimitive
{}

/// Collection of methods needed to implement maths in
/// GF(2<sup>x</sup>).  Provides (slow) default implementations that
/// can be overriden by optimised versions.
pub trait GaloisField {
    /// Natural storage class (u8, u16, u32, etc.) for storing
    /// elements of the field.
    type E  : ElementStore;
    
    /// The next largest natural storage class, eg if E is u8, then EE
    /// should be U16. Used for:
    ///
    /// * storing/returning field polynomial, which is always larger
    /// than the largest field element
    ///
    /// * storing the result of a non-modular (overflowing) multiply
    /// of two field elements
    type EE : ElementStore; // where Self::EE : From<Self::E>;
    /// As EE, but signed
    type SEE : ElementStore; // where Self::EE : From<Self::E>;

    /// Size of the field in bits, eg GF(2<sup>8</sup>) &rarr; 8
    const ORDER      : u16;
    /// High bit of field values, eg GF(2<sup>8</sup>) &rarr; 0x80
    const HIGH_BIT   : Self::E;
    /// High bit of field polynomial, eg GF(2<sup>8</sup>) &rarr; 0x100
    const POLY_BIT   : Self::EE;
    /// Mask for selecting all bits of field value,
    /// eg GF(2<sup>8</sup>) &rarr; 0xff
    const FIELD_MASK : Self::E;


    // If we try to implement one of the two *poly() methods in terms
    // of the other, we run into the problem of needing to convert
    // between E and EE, which is not worth the effort/hassle.
    
    /// The field polynomial with (implicit) high bit stripped off.
    fn poly(&self) -> Self::E;

    /// Field polynomial in full, with high poly bit set
    fn full_poly(&self) -> Self::EE;

    // Arithmetic operations
    
    /// Addition in Galois Fields GF(2<sup>x</sup>) is simply XOR
    fn add(&self, a : Self::E, b : Self::E) -> Self::E
    {
	a ^ b
    }
    /// Subtraction is the same as addition in GF(2<sup>x</sup>)
    fn sub(&self, a : Self::E, b : Self::E) -> Self::E
    {
	self.add(a, b)
    }

    /// Polynomial multiplication modulo the field polynomial
    fn mul(&self, mut a : Self::E, b : Self::E) -> Self::E
    {
	let poly = self.poly();
	//   let zero : Self::E = num::NumCast::from(0).unwrap();
	let zero = Self::E::zero();
	let one  = Self::E::one();
	//     let one  = Self::E::one();
	let field_mask : Self::E    = Self::FIELD_MASK;
	let high_bit : Self::E      = Self::HIGH_BIT;
	let mut result : Self::E    = if b & one != zero {a} else { zero };
	let mut bit : Self::E       = one + one;

	// mask needed in loop only to ensure GF(2**4) works correctly
	loop {
	    if a & high_bit != zero {
	 	// need to apply mask for GF(16)
	 	a = ((a << 1) ^ poly) & field_mask;
	    } else {
		// no mask here, as it can't overflow
		a = a << 1
	    }
	    if b & bit != zero {
	 	result = result ^ a
	    }
	    bit = bit << 1;
	    // 	// eprintln!("a: {}, b: {}, bit: {}, result: {}",
	    // 	//	  a,b,bit,result);
	    if bit & field_mask == zero { return result }
	}
	// unreachable!()
    }

    /// Division modulo to the field polynomial (default
    /// implementation calculated as a・b<sup>-1</sup>)
    fn div(&self, a : Self::E, b : Self::E) -> Self::E
    {
	self.mul(a, self.inv(b))
    }

    /// Calculate polynomial a<sup>-1</sup> modulo the field
    /// polynomial (using Extended Euclidean method)
    fn inv(&self, a : Self::E) -> Self::E
    {
	let (zero, one)    = (Self::E::zero(),
			      Self::E::one());
	let mut u : Self::E = self.poly();
	let mut v : Self::E = a;
	let mut z : Self::E = zero;
	let mut g : Self::E = one;
	let mut t : Self::E;	// always set before being used
	let mask = self.field_mask();
	// special cases of 1/1 and 1/0
	if a == zero || a == one { return a }
	// unroll first loop iteration (knowing initial i >= 0)
	// rustc can't determine that i is always positive here,
	// though, so we have to try_into()
	let fixup = if self.order() == 4 { 4 } else { 0 };
	let mut i : u32 = 1 + v.leading_zeros() - fixup;
	u = (u ^ v << i as usize) & mask;
	z = (z ^ g << i as usize) & mask;
	while u != one {
	    // in C: i=size_in_bits(u) - size_in_bits(v)
	    // size_in_bits is order - leading zeroes
	    // so we get i = (order - clz(u)) - (order - clz(v))
	    // the order cancels and we get:
	    // i = clz(v) - clz(u)
	    // further changed this to make i unsigned always
	    // i = v.leading_zeros() - u.leading_zeros();
	    if u.leading_zeros() > v.leading_zeros() {
		t=u; u=v; v=t;
		t=z; z=g; g=t;
		//		    i = u.leading_zeros() - v.leading_zeros();
	    } else {
		//		    i = v.leading_zeros() - u.leading_zeros();
	    }
	    i = v.leading_zeros() - u.leading_zeros();
	    u = (u ^ v << i as usize) & mask;
	    z = (z ^ g << i as usize) & mask;
	}
	z
    }
    
    /// Calculate polynomial a<sup>b</sup> modulo the field
    /// polynomial
    fn pow(&self, a : Self::E, mut b : Self::EE) -> Self::E
    where Self::E : Into<Self::EE> {
	let mut result : Self::E = a;
	let zero = Self::EE::zero();
	let one  = Self::E::one();
	let mut mask : Self::EE;

	// Enabling b % (field_size - 1). I'm doing this through
	// repeated subtraction rather than engage in gymnastics
	// involving type conversion from E -> usize -> E


	while b >=  (Self::FIELD_MASK).into() {
	    b = b - (Self::FIELD_MASK).into()
	}

	// identity below only works on field
	if b == zero
	// || b == self.field_mask().into()
	{ return one }
	// shift mask right if there are leading zeros
	// fixup for GF(16) to ignore unused top nibble

	// Since we have mask : EE, subtract any extra leading 0

	// want to move mask right so that it coincides with
	// highest bit of b. This alternative way works:
	let clz  = b.leading_zeros() as usize;
	let bits = zero.leading_zeros() as usize;
	//
	mask = Self::EE::one() << (bits - clz - 1);

	loop {
	    mask = mask >> 1;
	    if mask == zero { break }
	    result = self.mul(result, result);
	    if b & mask != zero { result = self.mul(result, a) }
	}
	result
    }

    // Long multiplication and modular reduction are not as useful to
    // users, though it is handy to have a default implementation for
    // adapted implementations. I can't make these private, but I have
    // options:
    //
    // * preface method name with _ as a *hint* that this is not
    //   usually not meant to be called by a user
    //
    // * hide the function in the user documentation

    #[doc(hidden)]
    /// Long Polynomial multiplication, *not* modulo the field polynomial

    // Note the extra type constraint for conversion E -> EE
    // (I could do something similar for poly/full_poly)
    fn mull(&self, a : Self::E, b : Self::E) -> Self::EE
    where Self::EE: From<Self::E>
    {
	let ezero  = Self::E::zero();
	let eezero = Self::EE::zero();

	// shift a, mask against b
	let mut aa : Self::EE = a.into();
	let mut res           = eezero;
	let mut mask          = Self::E::one();
	loop {
	    if b & mask != ezero { res = res ^ aa }
	    mask = mask << 1;
	    aa   = aa   << 1;
	    if aa == eezero      { return res }
	}
    }

    #[doc(hidden)]
    /// Bitwise modular reduction from EE to E

    // Here we have the reverse problem of converting an EE to E
    fn mod_reduce(&self, mut a : Self::EE) -> Self::E
    where Self::E: From<Self::EE>
    {
	let eezero   = Self::EE::zero();
	let mut poly = self.full_poly()  << (Self::ORDER - 1).into();
	let mut mask = Self::POLY_BIT    << (Self::ORDER - 1).into();
        loop {
	    if a & mask != eezero  { a = a ^ poly    }
	    if a < Self::POLY_BIT  { return a.into() }
	    mask = mask >> 1;
	    poly = poly >> 1;
	}
    }

    // Vector operations
    fn vec_sum_elements(&self, v : &[Self::E]) -> Self::E {
	let mut sum = Self::E::zero();
	for e in v.iter() {
            sum = sum ^ *e;
	}
	sum
    }

    fn vec_add_vec_in_place(&self,
			    dest  : &mut [Self::E],
			    other : &[Self::E] ) {
	assert_eq!(dest.len(), other.len());
	for (d,o) in dest.iter_mut().zip(other) {
	    *d = *d ^ *o
	}
    }

    fn vec_add_vecs_giving_other(&self,
				 dest  : &mut [Self::E],
				 a : &[Self::E],
				 b : &[Self::E]) {
	assert_eq!(dest.len(), a.len());
	assert_eq!(a.len(), b.len());
	let (mut a_iter, mut b_iter) = (a.iter(), b.iter());
	for d in dest.iter_mut() {
	    *d = *a_iter.next().unwrap() ^ *b_iter.next().unwrap()
	}
    }

    fn vec_cross_product(&self,
			    dest : &mut [Self::E],
			    a : &[Self::E],
			    b : &[Self::E] ) {
	assert_eq!(dest.len(), a.len());
	assert_eq!(a.len(), b.len());

	let (mut a_iter, mut b_iter) = (a.iter(), b.iter());
	for d in dest.iter_mut() {
            *d = self.mul(*a_iter.next().unwrap(), *b_iter.next().unwrap())
	}
    }

    fn vec_dot_product(&self, a : &[Self::E], b : &[Self::E]) -> Self::E
    {
	assert_eq!(a.len(), b.len());
	let mut sum = Self::E::zero();
	for (a_item, b_item) in a.iter().zip(b) {
            sum = sum ^ self.mul(*a_item, *b_item);
	}
	sum
    }

    fn vec_constant_scale_in_place(&self,
				   dest : &mut [Self::E],
				   a    : Self::E) {
	for d in dest.iter_mut() {
            *d = self.mul(*d,  a)
	}
    }

    // Obviously, this should take three *slices* ...
    fn vec_fma_in_place(&self,
			dest : &mut [Self::E],
			a    : Self::E,
			b    : Self::E ) {
	for d in dest.iter_mut() {
            *d = self.mul(*d, a) ^ b;
	}
    }


    // Other accessors provide syntactic sugar
    
    /// Access Self::HIGH_BIT as a method
    fn high_bit(&self)  -> Self::E { Self::HIGH_BIT   }

    /// Access Self::ORDER as a method
    fn order(&self)     -> u16     { Self::ORDER      }

    /// Access Self::FIELD_MASK as a method
    fn field_mask(&self)-> Self::E { Self::FIELD_MASK }


}


/// A type implementing (default) maths in GF(2<sup>4</sup>)
pub struct F4  { pub full : u8,  pub compact : u8 }

/// A type implementing (default) maths in GF(2<sup>8</sup>)
pub struct F8  { pub full : u16, pub compact : u8 }
	
/// A type implementing (default) maths in GF(2<sup>16</sup>)
pub struct F16 { pub full : u32, pub compact : u16 }

/// A type implementing (default) maths in GF(2<sup>32</sup>)
pub struct F32 { pub full : u64, pub compact : u32 }

impl GaloisField for F4 {
    type E = u8;
    type EE = u8;
    type SEE = i8;

    // we have to redeclare types for constants
    const ORDER      : u16 = 4;
    const POLY_BIT   : u8  = 0x10;
    const FIELD_MASK : u8  = 0x0f;
    const HIGH_BIT   : u8  = 0x08;

    // the two required methods (everything else is default)
    fn poly(&self)      -> u8  { self.compact }
    fn full_poly(&self) -> u8  { self.full }
}

// Constructor for GF(2<sup>4</sup>)
#[allow(dead_code)]
/// Create a new GF(2<sup>4</sup>) field with a supplied field
/// polynomial (using the default implementation)
pub fn new_gf4(full : u8, compact : u8) -> F4  {
    F4 {full, compact}
}

impl GaloisField for F8 {
    type E = u8;
    type EE = u16;
    type SEE = i16;

    // we have to redeclare types for constants
    const ORDER      : u16 = 8;
    const POLY_BIT   : u16 = 0x100;
    const FIELD_MASK : u8  = 0xff;
    const HIGH_BIT   : u8  = 0x80;

    // the two required methods (everything else is default)
    fn poly(&self)      -> u8  { self.compact }
    fn full_poly(&self) -> u16 { self.full }
}

// Constructor for GF(2<sup>8</sup>)
#[allow(dead_code)]
/// Create a new GF(2<sup>8</sup>) field with a supplied field
/// polynomial (using the default implementation)
pub fn new_gf8(full : u16, compact : u8) -> F8  {
    F8 {full, compact}
}

impl GaloisField for F16 {
    type E = u16;
    type EE = u32;
    type SEE = i32;

    // we have to redeclare types for constants
    const ORDER      : u16 = 8;
    const POLY_BIT   : u32 = 0x10000;
    const FIELD_MASK : u16 = 0xffff;
    const HIGH_BIT   : u16 = 0x8000;

    // the two required methods (everything else is default)
    fn poly(&self)      -> u16 { self.compact }
    fn full_poly(&self) -> u32 { self.full }
}

// Constructor for GF(2<sup>16</sup>)
#[allow(dead_code)]
/// Create a new GF(2<sup>16</sup>) field with a supplied field
/// polynomial (using the default implementation)
pub fn new_gf16(full : u32, compact : u16) -> F16  {
    F16 { full, compact }
}

impl GaloisField for F32 {
    type E = u32;
    type EE = u64;
    type SEE = i64;

    // we have to redeclare types for constants
    const ORDER      : u16 = 32;
    const POLY_BIT   : u64 = 0x100000000;
    const FIELD_MASK : u32 = 0xffffffff;
    const HIGH_BIT   : u32 = 0x80000000;

    // the two required methods (everything else is default)
    fn poly(&self)      -> u32 { self.compact }
    fn full_poly(&self) -> u64 { self.full }
}

// Constructor for GF(2<sup>32</sup>)
#[allow(dead_code)]
/// Create a new GF(2<sup>32</sup>) field with a supplied field
/// polynomial (using the default implementation)
pub fn new_gf32(full : u64, compact : u32) -> F32  {
    F32 { full, compact }
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn new_gf4_works() {
	let f = new_gf4(19, 3);
	let t = f.mul(1,1);
	assert_eq!(t, 1)
    }

    #[test]
    fn new_gf8_works() {
	let f = new_gf8(0x11d, 0x1d);
	let t = f.mul(1,1);
	assert_eq!(t, 1)
    }

    #[test]
    fn new_gf16_works() {
	let f = new_gf16(0x1_002b, 0x2b);
	let t = f.mul(1,1);
	assert_eq!(t, 1)
    }

    #[test]
    fn new_gf32_works() {
	let f = new_gf32(0x1_0000_008d, 0x8d);
	let t = f.mul(1,1);
	assert_eq!(t, 1)
    }

    #[test]
    fn zero_a_mod() {		// make "a" operand zero
	let obj = new_gf4(19, 3);
	let mut failed = 0;
	for i in 0..16 {
	    if obj.mul(0, i) != 0 { failed += 1 }
	}
	if failed > 0 { panic!("GF(16): failed {}/16 tests of a=0 * b", failed) }
    }

    #[test]
    fn zero_b_mod() {		// make "b" operand zero
	let obj = new_gf4(19, 3);
	let mut failed = 0;
	for i in 0..16 {
	    if obj.mul(i, 0) != 0 { failed += 1 }
	}
	if failed > 0 { panic!("GF(16): failed {}/16 tests of a * b=0", failed) }
    }

    #[test]
    fn one_a_mod() {		// make "a" operand one
	let obj = new_gf4(19, 3);
	let mut failed = 0;
	for i in 0..16 {
	    if obj.mul(i, 1) != i { failed += 1 }
	}
	if failed > 0 { panic!("GF(16): failed {}/16 tests of a=1 * b", failed) }
    }

    #[test]
    fn one_b_mod() {		// make "b" operand one
	let obj = new_gf4(19, 3);
	let mut failed = 0;
	for i in 0..16 {
	    if obj.mul(1, i) != i { failed += 1 }
	}
	if failed > 0 { panic!("GF(16): failed {}/16 tests of a * b=1", failed) }
    }

    // test that ab = ba for all a, b in {2..15}
    #[test]
    fn commutative_mod1() {
	let mut failed = 0;
	let obj = new_gf4(19, 3);
	for i in 2..16 {
	    for j in 2..16 {
		let result = obj.mul(i,j);
		if result != obj.mul(j, i) { failed += 1 }
	    }
	}
	if failed > 0 { panic!("GF(16): failed {}/169 commutativity tests", failed) }
    }

    // We can actually test "straight" multiply so long as we pick
    // operands that don't overflow, like these
    #[test]
    fn x_times_x1_mod1_straight() {
	let obj = new_gf4(19, 3);
	let (a, b) = (0b0000_0010, 0b0000_0011);
	assert_eq!(obj.mul(a,b), 0b0000_0110, "straight x(x+1)");
    }
    
    // another non-overflowing case
    #[test]
    fn xx_times_x1_mod1_straight() {
	let obj = new_gf4(19, 3);
	let (a, b) = (0b0000_0100, 0b0000_0011);
	assert_eq!(obj.mul(a,b), 0b0000_1100, "straight xx(x+1)");
    }

    #[test]
    fn distributive_1() {
	let obj = new_gf4(19, 3);
	let (a, b) = (0b0000_0011, 0b0000_1010);
	// first just call multiply on a, b
	let result = obj.mul(a, b);

	// now check that the distributive property held
	let check =  obj.add(obj.mul(0b0000_0010, b),
			     obj.mul(0b0000_0001, b));
	assert_eq!(result, check, "distributive 1");
    }

    #[test]
    fn distributive_2() {
	let obj = new_gf4(19, 3);
	let (a, b) = (0b0000_1010, 0b0000_1111);
	// first just call multiply on a, b
	let result = obj.mul(a, b);

	// now check that the distributive property held
	let check =  obj.add(obj.mul(0b0000_1000, b),
			     obj.mul(0b0000_0010, b));
	assert_eq!(result, check, "distributive 2");
    }

    // This test is about generators, but it only calls pow()
    #[test]
    fn generator_should_loop () {
	// values emitted from known_fields_have_generators() test
	let poly_19_known_gs = vec![2u8, 3, 4, 5, 9, 11, 13, 14];
	let poly_25_known_gs = vec![2u8, 4, 6, 7, 9, 12, 13, 14];
	let poly_31_known_gs = vec![3u8, 5, 6, 7, 9, 10, 11, 14];

	// writing a new bootstrap_mod_power fn for this, which I will
	// test separately. For now, make it a stub that causes this
	// test to pass!

	for (poly, gens) in [(19u8, poly_19_known_gs),
			     (25u8, poly_25_known_gs),
			     (31u8, poly_31_known_gs) ].iter() {
	    // if the generator "wraps around" correctly then:
	    // g ** 15 = 1
	    // g ** 16 = g ** 1 = g

	    for g in gens {
		let obj = new_gf4(*poly, *poly - 16);
		assert_eq!(1, obj.pow(*g, 15),
			   "poly,generator ({}, {}) failed g ** 15 = 1", *poly, *g);
		// can't call second part of test because it messes up
		// clz() call in pow for GF(16). Actually, this is a
		// problem in general, since I pass second operand as
		// type T, when it should really be P. Need a more
		// robust clz() that is aware of word size and field
		// size.  I don't feel like fixing this now because
		// it's going to open up the type conversion can of
		// worms again.
		
		assert_eq!(*g, obj.pow(*g, 16),
			   "poly,generator ({}, {}) failed g ** 16 = g", *poly, *g);
	    }
	}
    }


    #[test]
    fn mod_power_15() { // was failing
	let obj = new_gf4(19, 3);
	// would fail for a=0 because of corner case:
	// 0**0 = 0**15 = 1
	// (but 0**5 = 0!)
	for a in 1..=15 {
	    let direct = obj.pow(a, 15);

	    let indirect_5 = obj.pow(a, 5);
	    let indirect_15 = obj.pow(indirect_5, 3);

	    assert_eq!(direct, indirect_15, "a = {}", a);
	}
    }
    
    #[test]
    fn mod_power_16 () { // was failing, but note that we still can't
	// raise to 16th power directly. Code changes needed to
	// support that.

	// Next, use the decomposition 16 = 2**4 = 4**2
	let obj = new_gf4(19, 3);
	for a in 0..=15 {
	    let a2nd = obj.pow(a, 2);
	    let a4th = obj.pow(a, 4);

	    assert_eq!(a4th, obj.pow(a2nd, 2));

	    let ans_1 = obj.pow(a4th, 4);

	    // other decomposition a**4 * a**4 * a**4 * a**4
	    let a8th = obj.mul(a4th, a4th);
	    let a16th = obj.pow(a8th, 2);

	    assert_eq!(ans_1, a16th);
	}
    }

    #[test]
    fn mod_power_16_directly () {
	// Originally, couldn't raise a field to a higher power than
	// could be stored as a field element. Now that is possible,
	// so test it.

	let obj = new_gf4(19, 3);
	for a in 2..=15 {
	    let a15th = obj.pow(a, 15);
	    let a16th = obj.pow(a, 16); // was impossible before

	    assert_eq!(a16th, obj.mul(a, a15th), "fail for pow({},15/16)", a);
	}
    }

    // #[test]
    // see why power isn't working ...
    fn _debug_pow() {
	eprintln!("Testing 2**b for b 0..=15");
	let obj = new_gf4(19, 3);
	for i in 0..=16 {
	    let res = obj.pow(2, i);
	    eprintln!("2**{} = {}", i, res)
	}
    }

    
    // For order 8, I can copy/paste some tests from gf_2p8 crate
    #[test]
    fn test_new_gf8() {
	let obj = F8 { full : 0x11b, compact : 0x1b };
	let test = obj.mul(1,1);
	assert_eq!(test, 1);
    }

    #[test]
    fn new_field_from_new() {
	let f = new_gf8(0x11b, 0x1b);
	assert_eq!(8, F8::ORDER, "access ORDER associated type");
	assert_eq!(0x11b, f.full_poly(), "access full_poly method");
    }    
    
    #[test]
    fn test_multiply_53_ca() {
	assert_eq!(1, new_gf8(0x11b, 0x1b).mul(0x53, 0xCA))
    }
    
    #[test]
    fn test_multiply_00_00() {
	assert_eq!(0, new_gf8(0x11b, 0x1b).mul(0x00, 0x00))
    }
    
    #[test]
    fn test_multiply_00_01() {
	assert_eq!(0, new_gf8(0x11b, 0x1b).mul(0x00, 0x01))
    }
    
    #[test]
    fn test_multiply_01_01() {
	assert_eq!(1, new_gf8(0x11b, 0x1b).mul(0x01, 0x01))
    }
    
    #[test]
    fn test_multiply_01_02() {
	assert_eq!(2, new_gf8(0x11b, 0x1b).mul(0x01, 0x02))
    }

    #[test]
    fn test_inv_01() {
	assert_eq!(1, new_gf8(0x11b, 0x1b).inv(0x01))
    }

    #[test]
    fn test_inv_53() {
	assert_eq!(0xca, new_gf8(0x11b, 0x1b).inv(0x53))
    }

    #[test]
    fn test_inv_ca() {
	assert_eq!(0x53, new_gf8(0x11b, 0x1b).inv(0xca))
    }

    #[test]
    fn test_div_01_53() {
	assert_eq!(0xca, new_gf8(0x11b, 0x1b).div(1,0x53))
    }

    #[test]
    fn test_div_01_ca() {
	assert_eq!(0x53, new_gf8(0x11b, 0x1b).div(1,0xca))
    }

    // for powers, we (I) know that 3 is a generator, so:
    // 3**0 == 3**255 = 1
    // 3**-1 == 3**254 = inv(3)
    //
    // we can't pass in negative powers, though...
    #[test]
    fn test_pow_3_0() {
	assert_eq!(0x01, new_gf8(0x11b, 0x1b).pow(3,0))
    }
    
    #[test]
    fn test_pow_3_255() {
	assert_eq!(0x1, new_gf8(0x11b, 0x1b).pow(3,255))
    }
    
    #[test]
    fn test_pow_3_254() {
	let f = new_gf8(0x11b, 0x1b);
	assert_eq!(f.pow(3,254), f.inv(3))
    }

    // compare power with multiplication
    #[test]
    fn test_pow_7_2() {
	let f = new_gf8(0x11b, 0x1b);
	assert_eq!(f.pow(7,2), f.mul(7,7))
    }
    
    #[test]
    fn test_pow_7_3() {
	let f = new_gf8(0x11b, 0x1b);
	assert_eq!(f.pow(7,3), f.mul(f.pow(7,2),7))
    }
    
    // special case: 0**0 = 1
    #[test]
    fn test_pow_0_0() {
	let f = new_gf8(0x11b, 0x1b);
	assert_eq!(1,f.pow(0,0))
    }

    // special case: 1/0 = 0 (undefined behaviour!)
    #[test]
    fn test_div_0_0() {
	let f = new_gf8(0x11b, 0x1b);
	assert_eq!(0,f.div(0,0))
    }

    // related cases: non-zero / zero = 0 (also UB)
    #[test]
    fn test_div_i_0() {
	let f = new_gf8(0x11b, 0x1b);
	for i in 1..=255 {
	    assert_eq!(0,f.div(i,0))
	}
    }

    // additional: make sure pow(a,256), pow(a,257) works in this
    // field
    #[test]
    fn test_pow_257() {
	let f = F8 { full : 0x11b, compact : 0x1b };
	let a = f.pow(12,254);
	let b = f.pow(12,255);
	let c = f.pow(12,256);
	let d = f.pow(12,257);
	assert_eq!(b,f.mul(a,12));
	assert_eq!(c,f.mul(b,12));
	assert_eq!(d,f.mul(c,12));
    }

    #[test]
    fn test_new_gf16() {
	let obj = new_gf16(0x1002b, 0x2b);
	let test = obj.mul(1,1);
	assert_eq!(test, 1);
    }

    #[test]
    fn test_new_gf32() {
	let obj = new_gf32(0x10000008d, 0x8d);
	let test = obj.mul(1,1);
	assert_eq!(test, 1);
    }

    #[test]
    fn access_assoc_type() {
	type F = <F8 as GaloisField>::E;
	let _a : F = 1;
	assert_eq!(1, std::mem::size_of::<F>());
    }

    // getting access to the internal associated type is
    // not easy. Giving up for now.
    //
    // "cannot provide explicit generic arguments when `impl Trait` is
    // used in argument position"

    fn test_impl_assoc_type<T>(_f : & T, s : usize)
    where T : GaloisField
    {
	//type Y = T::FieldStore;
	let one : T::E = T::E::one();
	let mut _a : T::E = num::NumCast::from(1).unwrap();
	// cannot infer type:
	// let mut a : T::FieldStore = From::from(1);
	_a = _a + one;
	assert_eq!(s, std::mem::size_of::<T::E>());
    }

    #[test]
    fn call_impl_assoc_type() {
	let f4 = new_gf4(29, 13);
	test_impl_assoc_type(&f4, 1);
	let f8 = new_gf4(29, 13);
	test_impl_assoc_type(&f8, 1);
	let f16 = new_gf16(29, 13);
	test_impl_assoc_type(&f16, 2);
	let f32 = new_gf32(29, 13);
	test_impl_assoc_type(&f32, 4);
    }

    // Make sure that long multiply + mod reduce agree with regular
    // mul method
    #[test]
    fn long_mul_mod_reduce_conformance() {
	let f = new_gf4(19, 3);
	for a in 0..=15 {
	    for b in 0..=15 {
		let longmul = f.mull(a,b);
		assert_eq!(f.mul(a,b), f.mod_reduce(longmul));
	    }
	}
    }

    #[test]
    fn test_gf8_inverses() {
	let f = new_gf8(0x11b, 0x1b);
	for a in 0..=255 {
	    let inv = f.inv(a);
	    let invinv = f.inv(inv);
	    assert_eq!(a, invinv);
	}
    }

    #[test]
    fn access_lmull() {
	use crate::tables::mull;
	assert_eq!(mull::RMULL.len(), 4096);
    }
}
