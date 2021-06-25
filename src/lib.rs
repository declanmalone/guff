//! # Grand Unified Finite Field library
//!
//! Implements GF(2<sup>x</sup>) for various natural sizes (currently
//! 8, 16 and 32).
//! 
//! 
//! 


use num_traits;
use num::{PrimInt,One,Zero};


pub trait FieldStore : 'static + Copy
    + num_traits::int::PrimInt
    + std::fmt::Display
    + num::FromPrimitive + num::ToPrimitive
    // + num::FromPrimitive + num::ToPrimitive
    // + num::Integer
    // + num::cast::AsPrimitive<Self::E>
    // + num::Integer + num::One
{}
impl<T> FieldStore for T
where T : 'static + Copy
    + num_traits::int::PrimInt
    + std::fmt::Display
    + num::FromPrimitive + num::ToPrimitive
{}

pub trait GaloisField {
    /// Natural storage class (u8, u16, u32, etc.) for storing
    /// elements of the field.
    type E  : FieldStore;
    /// The next largest natural storage class, eg if E is u8, then EE
    /// should be U16. Used for:
    ///
    ///* storing/returning field polynomial, which is always larger
    ///  than the largest field element
    ///
    ///* storing the result of a non-modular (overflowing) multiply of
    /// two field elements
    type EE : FieldStore;

    /// The field polynomial with (implicit) high bit stripped off.
    fn poly(&self) -> Self::E;
    /// Field polynomial in full
    fn full_poly(&self) -> Self::EE;
    /// These will become associated constants
    fn high_bit() -> Self::E;
    fn order() -> u16;
    fn field_mask() -> Self::E;

    fn add(&self, a : Self::E, b : Self::E) -> Self::E
    {
	a ^ b
    }

    // upgrade storage class to prevent overflow
    // fn mul_no_mod(&self, a : Self::E, b : Self::E)
    // 	      -> Self::EE
    // {
    //     // implement straight, non-modular addition
    //     6u8.into()
    // }
    fn mul(&self, mut a : Self::E, b : Self::E) -> Self::E
    {
	let poly = self.poly();
	//   let zero : Self::E = num::NumCast::from(0).unwrap();
	let zero : Self::E = num::zero();
	let one  : Self::E = num::one();
	//     let one  = Self::E::one();
	let field_mask : Self::E    = Self::field_mask();
	let high_bit : Self::E      = Self::high_bit();
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

    fn div(&self, a : Self::E, b : Self::E) -> Self::E
    {
	self.mul(a, self.inv(b))
    }
    
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
	let one : Self::E        = num::one();
	let mut mask : Self::EE;

	// identity below only works on field
	if b == num::zero()
	// || b == self.field_mask()
	{ return one }
	// shift mask right if there are leading zeros
	// fixup for GF(16) to ignore unused top nibble
	
	// Since we have mask : P, subtract any extra leading 0

	// want to move mask right so that it coincides with
	// highest bit of b. This alternative way works:
	let clz = b.leading_zeros() as usize;
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
    
}


struct F8;








#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
