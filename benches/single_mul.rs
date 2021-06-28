

use guff::{GaloisField, new_gf4, new_gf8};
use guff::good::{new_gf4_0x13, new_gf8_0x11b};


// const ref_f : F4 = F4 {full : 19, compact : 3};

// fn single_mul(a : u8, b : u8)
// {
//     ref_f.mul(a,b);
// }



// despite docs, should have no main() here.
// #![allow(unused)]
//fn main() {

//    

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use criterion::BenchmarkId;

// Use like-for-like harness (bench_with_input instead of
// bench_function)
//
// Also, multiply all pairs of values

pub fn ref_gf4_mul(c: &mut Criterion) {
    let ref_f = new_gf4(19,3);
    c.bench_with_input(
		       BenchmarkId::new("gf4 mul", "ref"),
		       &ref_f,
		       |b, f| {
			   b.iter(||
				  for i in 0..=15 {
				      for j in 0..=15 {
					  f.mul(i,j);
				      }
				  }
			   );
		       });
}
    
pub fn ref_gf4_mull(c: &mut Criterion) {
    let ref_f = new_gf4(19,3);
    c.bench_with_input(
		       BenchmarkId::new("gf4 mul", "mull"),
		       &ref_f,
		       |b, f| {
			   b.iter(||
				  for i in 0..=15 {
				      for j in 0..=15 {
					  f.mull(black_box(i),j);
				      }
				  }
			   );
		       });
}
    
pub fn ref_gf4_mull_reduce(c: &mut Criterion) {
    let ref_f = new_gf4(19,3);
    c.bench_with_input(
		       BenchmarkId::new("gf4 mul", "mull-reduce"),
		       &ref_f,
		       |b, f| {
			   b.iter(||
				  for i in 0..=15 {
				      for j in 0..=15 {
					  f.mod_reduce(f.mull(i,j));
				      }
				  }
			   );
		       });
}
    
pub fn good_gf4_mul(c: &mut Criterion) {
    let good_f4 = new_gf4_0x13();
    c.bench_with_input(
		       BenchmarkId::new("gf4 mul", "good"),
		       &good_f4,
		       |b, f| {
			   b.iter(||
				  for i in 0..=15 {
				      for j in 0..=15 {
					  f.mul(i,j);
				      }
				  }
			   );
			   });
}

pub fn ref_gf4_inv(c: &mut Criterion) {
    let ref_f = new_gf4(19,3);
    c.bench_with_input(
		       BenchmarkId::new("gf4 inv", "ref"),
		       &ref_f,
		       |b, f| {
			   b.iter(||
				  for i in 0..=15 {
					  f.inv(i);
				  }
			   );
		       });
}

pub fn good_gf4_inv(c: &mut Criterion) {
    let good_f = new_gf4_0x13();
    c.bench_with_input(
		       BenchmarkId::new("gf4 inv", "good"),
		       &good_f,
		       |b, f| {
			   b.iter(||
				  for i in 0..=15 {
					  f.inv(i);
				  }
			   );
		       });
}

// gf(256)
pub fn ref_gf8_mul(c: &mut Criterion) {
    let ref_f = new_gf8(0x11b,0x1b);
    c.bench_with_input(
		       BenchmarkId::new("gf8 mul", "ref"),
		       &ref_f,
		       |b, f| {
			   b.iter(||
				  for i in 0..=255 {
				      for j in 0..=255 {
					  f.mul(i,j);
				      }
				  }
			   );
		       });
}

pub fn good_gf8_mul(c: &mut Criterion) {
    let good_f8 = new_gf8_0x11b();
    c.bench_with_input(
		       BenchmarkId::new("gf8 mul", "good"),
		       &good_f8,
		       |b, f| {
			   b.iter(||
				  for i in 0..=255 {
				      for j in 0..=255 {
					  f.mul(i,j);
				      }
				  }
			   );
		       });
}



criterion_group!(benches,
		 ref_gf4_mul, good_gf4_mul,
		 ref_gf4_mull, ref_gf4_mull_reduce, 
		 ref_gf4_inv, good_gf4_inv,
		 ref_gf8_mul, good_gf8_mul
);
criterion_main!(benches);

//}
