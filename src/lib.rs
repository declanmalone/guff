//! # Grand Unified Finite Field library*
//!
//! Implements GF(2<sup>x</sup>) for various "natural" sizes
//! such as 2<sup>8</sup>, 2<sup>16</sup> and 2<sup>32</sup>.
//!
//! My goals for this crate are to:
//! 
//! 1. help me learn to write good modules in Rust
//! 
//! 2. help interested users learn about finite fields (ie, Galois
//! fields)
//! 
//! 3. provide a generic baseline implementation of basic maths
//! (add, multiply, divide, etc.) over finite fields
//! 
//! 4. explore various optimisations/adaptations (including
//! table-based lookups and architecture-specific SIMD code) that can
//! selectively override some/all of the default implementations
//! (while remaining compatible with other implementations).
//! 
//! Also to:
//! 
//! 5. provide some useful utility functions that go beyond just
//! `add`, `mul`, `div`, etc. (eg, determining whether a field
//! polynomial is primitive, or generating lookup tables for different
//! kinds of optimisations)
//! 
//! See the top-level `Readme` for more information about the above.
//! 
//! # Basic Use: doing maths in a particular field
//! 
//! As a user, the steps to take to use this library are:
//!
//! * decide what "class" of field you want to use (GF(2<sup>8</sup>),
//! GF(2<sup>16</sup>), etc.);
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








//! # Crate Name
//! 
//! \* The crate name is deliberately hyperbolic:
//!
//! > Noun *guff* - unacceptable behavior (especially ludicrously false statements)

use num_traits;
use num::{PrimInt,One,Zero};

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

/// 
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
	let zero : Self::E = num::zero();
	let one  : Self::E = num::one();
	//     let one  = Self::E::one();
	let field_mask : Self::E    = Self::FIELD_MASK;
	let high_bit : Self::E      = Self::HIGH_BIT;
	let mut result : Self::E    = if b & one != zero {a} else { zero };
	let mut bit : Self::E       = one + one;
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
    /// implementation calculated as a・b<sup>-1</sup>
    fn div(&self, a : Self::E, b : Self::E) -> Self::E
    {
	self.mul(a, self.inv(b))
    }

    /// Calculate polynomial a・b<sup>-1</sup> modulo the field
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
	// special cases of 1/1 and 1/0
	if a == zero || a == one { return a }
	// unroll first loop iteration (knowing initial i >= 0)
	// rustc can't determine that i is always positive here,
	// though, so we have to try_into()
	let mut i : u32 = 1 + v.leading_zeros();
	u = u ^ v << i as usize;
	z = z ^ g << i as usize;
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
	    u = u ^ v << i as usize;
	    z = z ^ g << i as usize;
	}
	z
    }
    
    fn pow(&self, a : Self::E, b : Self::EE) -> Self::E {
	let mut result : Self::E = a;
	let zero : Self::EE     = num::zero();
	let one  : Self::E      = num::one();
	let mut mask : Self::EE;

	// identity below only works on field
	if b == num::zero()
	// || b == self.field_mask()
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
    where <Self as GaloisField>::EE: From<<Self as GaloisField>::E>
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

    
    // Other accessors provide syntactic sugar
    
    /// Access Self::HIGH_BIT as a method
    fn high_bit(&self)  -> Self::E { Self::HIGH_BIT   }

    /// Access Self::ORDER as a method
    fn order(&self)     -> u16     { Self::ORDER      }

    /// Access Self::FIELD_MASK as a method
    fn field_mask(&self)-> Self::E { Self::FIELD_MASK }


}


// 
pub struct F4{pub poly : u8}
pub struct F8{pub poly : u8}









#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
