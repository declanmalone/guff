//! # Long, non-modular polynomial multiplication lookup table
//!
//! The functions here help when multiplying polynomials that are
//! larger than `u8` in size. They are backed by an 8k static lookup
//! table containing all possible polynomial multiplications between a
//! 4-bit polynomial fragment and an 8-bit polynomial fragment.
//!
//! Example: multiplying a 32-bit polynomial by an 8-bit one
//!
//! ``` rust
//!
//! use guff::tables::mull::{rmull};
//!
//! let a : u32 = 0x1234_5678;
//! let b : u8  = 0x9a;
//!
//! // break a into separate bytes, b into nibbles
//! let a3 : u8 = ((a >> 24) as u8);
//! let a2 : u8 = ((a >> 16) as u8);
//! let a1 : u8 = ((a >>  8) as u8);
//! let a0 : u8 = ((a      ) as u8);
//!
//! let b1 = b >> 4;
//! let b0 = b & 0x0f;
//!
//! let mut result : u64 = 0;
//!
//! // piece-wise multiplication using schoolbook method
//! result ^= (rmull(a0,b0) as u64) << 0  ;
//! result ^= (rmull(a1,b0) as u64) << 8  ;
//! result ^= (rmull(a2,b0) as u64) << 16 ;
//! result ^= (rmull(a3,b0) as u64) << 24 ;
//! result ^= (rmull(a0,b1) as u64) << 4  ;
//! result ^= (rmull(a1,b1) as u64) << 12 ;
//! result ^= (rmull(a2,b1) as u64) << 20 ;
//! result ^= (rmull(a3,b1) as u64) << 28 ;
//!
//! // compare result with GF(2**32) reference mull
//! use guff::{GaloisField,F32};
//!
//! assert_eq!(result,
//!            F32::mull(a.into(), b.into()));
//!
//!
//! ```

/// Straight multiply assuming small value is in low (rightmost) nibble position
#[inline(always)]
pub fn rmull(big : u8, small : u8) -> u16 {
    let index = ( ((small as u16) << 8) ^ (big as u16) ) as usize;
    let out = MULL[index];
    // eprintln!("Big {}, Small {}", big, small);
    // eprintln!("Index {}", index);
    // eprintln!("Out: {}", out);
    out
}

/// Straight multiply assuming small value is in high (leftmost) nibble position
#[inline(always)]
pub fn lmull(big : u8, small : u8) -> u16 {
    let index = ( ((small as u16) << 8) ^ (big as u16) ) as usize;
    let out = MULL[index];
    out << 4
}

/// Straight multiply with small value containing packed left, right nibbles

// Might be slightly slower than calling lmull, rmull directly, due to
// repeated calculation of shift/mask at the top, but the compiler
// should be able to optimise these out.
#[inline(always)]
pub fn lrmull(big : u8, small : u8) -> u16 {

    let l : u16 = (small >> 4   ).into();
    let r : u16 = (small & 0x0f ).into();

    // Exact same index calculations for l, r
    let l_index = ( (l << 8) ^ (big as u16) ) as usize;
    let l_out = MULL[l_index];

    let r_index = ( (r << 8) ^ (big as u16) ) as usize;
    let r_out = MULL[r_index];

    // Difference is that we left-shift the l result
    (l_out << 4) ^ r_out
}

// 12-bit multiply of nibble by 8-bit value
//
// RMULL = MULL[(byte << 4) | (nibble)] => byte * nibble
//
// I will reuse this table to calculate LMULL, which is the RMULL
// result shifted left 4 bits.

/// Lookup table for multiplying 8-bit poly fragment by 4-bit poly
/// fragment (straight multiplication; no modulus)
pub const MULL : [u16; 4096] = [
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 
  0, 1, 2, 3, 4, 5, 6, 7, 
  8, 9, 10, 11, 12, 13, 14, 15, 
  16, 17, 18, 19, 20, 21, 22, 23, 
  24, 25, 26, 27, 28, 29, 30, 31, 
  32, 33, 34, 35, 36, 37, 38, 39, 
  40, 41, 42, 43, 44, 45, 46, 47, 
  48, 49, 50, 51, 52, 53, 54, 55, 
  56, 57, 58, 59, 60, 61, 62, 63, 
  64, 65, 66, 67, 68, 69, 70, 71, 
  72, 73, 74, 75, 76, 77, 78, 79, 
  80, 81, 82, 83, 84, 85, 86, 87, 
  88, 89, 90, 91, 92, 93, 94, 95, 
  96, 97, 98, 99, 100, 101, 102, 103, 
  104, 105, 106, 107, 108, 109, 110, 111, 
  112, 113, 114, 115, 116, 117, 118, 119, 
  120, 121, 122, 123, 124, 125, 126, 127, 
  128, 129, 130, 131, 132, 133, 134, 135, 
  136, 137, 138, 139, 140, 141, 142, 143, 
  144, 145, 146, 147, 148, 149, 150, 151, 
  152, 153, 154, 155, 156, 157, 158, 159, 
  160, 161, 162, 163, 164, 165, 166, 167, 
  168, 169, 170, 171, 172, 173, 174, 175, 
  176, 177, 178, 179, 180, 181, 182, 183, 
  184, 185, 186, 187, 188, 189, 190, 191, 
  192, 193, 194, 195, 196, 197, 198, 199, 
  200, 201, 202, 203, 204, 205, 206, 207, 
  208, 209, 210, 211, 212, 213, 214, 215, 
  216, 217, 218, 219, 220, 221, 222, 223, 
  224, 225, 226, 227, 228, 229, 230, 231, 
  232, 233, 234, 235, 236, 237, 238, 239, 
  240, 241, 242, 243, 244, 245, 246, 247, 
  248, 249, 250, 251, 252, 253, 254, 255, 
  0, 2, 4, 6, 8, 10, 12, 14, 
  16, 18, 20, 22, 24, 26, 28, 30, 
  32, 34, 36, 38, 40, 42, 44, 46, 
  48, 50, 52, 54, 56, 58, 60, 62, 
  64, 66, 68, 70, 72, 74, 76, 78, 
  80, 82, 84, 86, 88, 90, 92, 94, 
  96, 98, 100, 102, 104, 106, 108, 110, 
  112, 114, 116, 118, 120, 122, 124, 126, 
  128, 130, 132, 134, 136, 138, 140, 142, 
  144, 146, 148, 150, 152, 154, 156, 158, 
  160, 162, 164, 166, 168, 170, 172, 174, 
  176, 178, 180, 182, 184, 186, 188, 190, 
  192, 194, 196, 198, 200, 202, 204, 206, 
  208, 210, 212, 214, 216, 218, 220, 222, 
  224, 226, 228, 230, 232, 234, 236, 238, 
  240, 242, 244, 246, 248, 250, 252, 254, 
  256, 258, 260, 262, 264, 266, 268, 270, 
  272, 274, 276, 278, 280, 282, 284, 286, 
  288, 290, 292, 294, 296, 298, 300, 302, 
  304, 306, 308, 310, 312, 314, 316, 318, 
  320, 322, 324, 326, 328, 330, 332, 334, 
  336, 338, 340, 342, 344, 346, 348, 350, 
  352, 354, 356, 358, 360, 362, 364, 366, 
  368, 370, 372, 374, 376, 378, 380, 382, 
  384, 386, 388, 390, 392, 394, 396, 398, 
  400, 402, 404, 406, 408, 410, 412, 414, 
  416, 418, 420, 422, 424, 426, 428, 430, 
  432, 434, 436, 438, 440, 442, 444, 446, 
  448, 450, 452, 454, 456, 458, 460, 462, 
  464, 466, 468, 470, 472, 474, 476, 478, 
  480, 482, 484, 486, 488, 490, 492, 494, 
  496, 498, 500, 502, 504, 506, 508, 510, 
  0, 3, 6, 5, 12, 15, 10, 9, 
  24, 27, 30, 29, 20, 23, 18, 17, 
  48, 51, 54, 53, 60, 63, 58, 57, 
  40, 43, 46, 45, 36, 39, 34, 33, 
  96, 99, 102, 101, 108, 111, 106, 105, 
  120, 123, 126, 125, 116, 119, 114, 113, 
  80, 83, 86, 85, 92, 95, 90, 89, 
  72, 75, 78, 77, 68, 71, 66, 65, 
  192, 195, 198, 197, 204, 207, 202, 201, 
  216, 219, 222, 221, 212, 215, 210, 209, 
  240, 243, 246, 245, 252, 255, 250, 249, 
  232, 235, 238, 237, 228, 231, 226, 225, 
  160, 163, 166, 165, 172, 175, 170, 169, 
  184, 187, 190, 189, 180, 183, 178, 177, 
  144, 147, 150, 149, 156, 159, 154, 153, 
  136, 139, 142, 141, 132, 135, 130, 129, 
  384, 387, 390, 389, 396, 399, 394, 393, 
  408, 411, 414, 413, 404, 407, 402, 401, 
  432, 435, 438, 437, 444, 447, 442, 441, 
  424, 427, 430, 429, 420, 423, 418, 417, 
  480, 483, 486, 485, 492, 495, 490, 489, 
  504, 507, 510, 509, 500, 503, 498, 497, 
  464, 467, 470, 469, 476, 479, 474, 473, 
  456, 459, 462, 461, 452, 455, 450, 449, 
  320, 323, 326, 325, 332, 335, 330, 329, 
  344, 347, 350, 349, 340, 343, 338, 337, 
  368, 371, 374, 373, 380, 383, 378, 377, 
  360, 363, 366, 365, 356, 359, 354, 353, 
  288, 291, 294, 293, 300, 303, 298, 297, 
  312, 315, 318, 317, 308, 311, 306, 305, 
  272, 275, 278, 277, 284, 287, 282, 281, 
  264, 267, 270, 269, 260, 263, 258, 257, 
  0, 4, 8, 12, 16, 20, 24, 28, 
  32, 36, 40, 44, 48, 52, 56, 60, 
  64, 68, 72, 76, 80, 84, 88, 92, 
  96, 100, 104, 108, 112, 116, 120, 124, 
  128, 132, 136, 140, 144, 148, 152, 156, 
  160, 164, 168, 172, 176, 180, 184, 188, 
  192, 196, 200, 204, 208, 212, 216, 220, 
  224, 228, 232, 236, 240, 244, 248, 252, 
  256, 260, 264, 268, 272, 276, 280, 284, 
  288, 292, 296, 300, 304, 308, 312, 316, 
  320, 324, 328, 332, 336, 340, 344, 348, 
  352, 356, 360, 364, 368, 372, 376, 380, 
  384, 388, 392, 396, 400, 404, 408, 412, 
  416, 420, 424, 428, 432, 436, 440, 444, 
  448, 452, 456, 460, 464, 468, 472, 476, 
  480, 484, 488, 492, 496, 500, 504, 508, 
  512, 516, 520, 524, 528, 532, 536, 540, 
  544, 548, 552, 556, 560, 564, 568, 572, 
  576, 580, 584, 588, 592, 596, 600, 604, 
  608, 612, 616, 620, 624, 628, 632, 636, 
  640, 644, 648, 652, 656, 660, 664, 668, 
  672, 676, 680, 684, 688, 692, 696, 700, 
  704, 708, 712, 716, 720, 724, 728, 732, 
  736, 740, 744, 748, 752, 756, 760, 764, 
  768, 772, 776, 780, 784, 788, 792, 796, 
  800, 804, 808, 812, 816, 820, 824, 828, 
  832, 836, 840, 844, 848, 852, 856, 860, 
  864, 868, 872, 876, 880, 884, 888, 892, 
  896, 900, 904, 908, 912, 916, 920, 924, 
  928, 932, 936, 940, 944, 948, 952, 956, 
  960, 964, 968, 972, 976, 980, 984, 988, 
  992, 996, 1000, 1004, 1008, 1012, 1016, 1020, 
  0, 5, 10, 15, 20, 17, 30, 27, 
  40, 45, 34, 39, 60, 57, 54, 51, 
  80, 85, 90, 95, 68, 65, 78, 75, 
  120, 125, 114, 119, 108, 105, 102, 99, 
  160, 165, 170, 175, 180, 177, 190, 187, 
  136, 141, 130, 135, 156, 153, 150, 147, 
  240, 245, 250, 255, 228, 225, 238, 235, 
  216, 221, 210, 215, 204, 201, 198, 195, 
  320, 325, 330, 335, 340, 337, 350, 347, 
  360, 365, 354, 359, 380, 377, 374, 371, 
  272, 277, 282, 287, 260, 257, 270, 267, 
  312, 317, 306, 311, 300, 297, 294, 291, 
  480, 485, 490, 495, 500, 497, 510, 507, 
  456, 461, 450, 455, 476, 473, 470, 467, 
  432, 437, 442, 447, 420, 417, 430, 427, 
  408, 413, 402, 407, 396, 393, 390, 387, 
  640, 645, 650, 655, 660, 657, 670, 667, 
  680, 685, 674, 679, 700, 697, 694, 691, 
  720, 725, 730, 735, 708, 705, 718, 715, 
  760, 765, 754, 759, 748, 745, 742, 739, 
  544, 549, 554, 559, 564, 561, 574, 571, 
  520, 525, 514, 519, 540, 537, 534, 531, 
  624, 629, 634, 639, 612, 609, 622, 619, 
  600, 605, 594, 599, 588, 585, 582, 579, 
  960, 965, 970, 975, 980, 977, 990, 987, 
  1000, 1005, 994, 999, 1020, 1017, 1014, 1011, 
  912, 917, 922, 927, 900, 897, 910, 907, 
  952, 957, 946, 951, 940, 937, 934, 931, 
  864, 869, 874, 879, 884, 881, 894, 891, 
  840, 845, 834, 839, 860, 857, 854, 851, 
  816, 821, 826, 831, 804, 801, 814, 811, 
  792, 797, 786, 791, 780, 777, 774, 771, 
  0, 6, 12, 10, 24, 30, 20, 18, 
  48, 54, 60, 58, 40, 46, 36, 34, 
  96, 102, 108, 106, 120, 126, 116, 114, 
  80, 86, 92, 90, 72, 78, 68, 66, 
  192, 198, 204, 202, 216, 222, 212, 210, 
  240, 246, 252, 250, 232, 238, 228, 226, 
  160, 166, 172, 170, 184, 190, 180, 178, 
  144, 150, 156, 154, 136, 142, 132, 130, 
  384, 390, 396, 394, 408, 414, 404, 402, 
  432, 438, 444, 442, 424, 430, 420, 418, 
  480, 486, 492, 490, 504, 510, 500, 498, 
  464, 470, 476, 474, 456, 462, 452, 450, 
  320, 326, 332, 330, 344, 350, 340, 338, 
  368, 374, 380, 378, 360, 366, 356, 354, 
  288, 294, 300, 298, 312, 318, 308, 306, 
  272, 278, 284, 282, 264, 270, 260, 258, 
  768, 774, 780, 778, 792, 798, 788, 786, 
  816, 822, 828, 826, 808, 814, 804, 802, 
  864, 870, 876, 874, 888, 894, 884, 882, 
  848, 854, 860, 858, 840, 846, 836, 834, 
  960, 966, 972, 970, 984, 990, 980, 978, 
  1008, 1014, 1020, 1018, 1000, 1006, 996, 994, 
  928, 934, 940, 938, 952, 958, 948, 946, 
  912, 918, 924, 922, 904, 910, 900, 898, 
  640, 646, 652, 650, 664, 670, 660, 658, 
  688, 694, 700, 698, 680, 686, 676, 674, 
  736, 742, 748, 746, 760, 766, 756, 754, 
  720, 726, 732, 730, 712, 718, 708, 706, 
  576, 582, 588, 586, 600, 606, 596, 594, 
  624, 630, 636, 634, 616, 622, 612, 610, 
  544, 550, 556, 554, 568, 574, 564, 562, 
  528, 534, 540, 538, 520, 526, 516, 514, 
  0, 7, 14, 9, 28, 27, 18, 21, 
  56, 63, 54, 49, 36, 35, 42, 45, 
  112, 119, 126, 121, 108, 107, 98, 101, 
  72, 79, 70, 65, 84, 83, 90, 93, 
  224, 231, 238, 233, 252, 251, 242, 245, 
  216, 223, 214, 209, 196, 195, 202, 205, 
  144, 151, 158, 153, 140, 139, 130, 133, 
  168, 175, 166, 161, 180, 179, 186, 189, 
  448, 455, 462, 457, 476, 475, 466, 469, 
  504, 511, 502, 497, 484, 483, 490, 493, 
  432, 439, 446, 441, 428, 427, 418, 421, 
  392, 399, 390, 385, 404, 403, 410, 413, 
  288, 295, 302, 297, 316, 315, 306, 309, 
  280, 287, 278, 273, 260, 259, 266, 269, 
  336, 343, 350, 345, 332, 331, 322, 325, 
  360, 367, 358, 353, 372, 371, 378, 381, 
  896, 903, 910, 905, 924, 923, 914, 917, 
  952, 959, 950, 945, 932, 931, 938, 941, 
  1008, 1015, 1022, 1017, 1004, 1003, 994, 997, 
  968, 975, 966, 961, 980, 979, 986, 989, 
  864, 871, 878, 873, 892, 891, 882, 885, 
  856, 863, 854, 849, 836, 835, 842, 845, 
  784, 791, 798, 793, 780, 779, 770, 773, 
  808, 815, 806, 801, 820, 819, 826, 829, 
  576, 583, 590, 585, 604, 603, 594, 597, 
  632, 639, 630, 625, 612, 611, 618, 621, 
  560, 567, 574, 569, 556, 555, 546, 549, 
  520, 527, 518, 513, 532, 531, 538, 541, 
  672, 679, 686, 681, 700, 699, 690, 693, 
  664, 671, 662, 657, 644, 643, 650, 653, 
  720, 727, 734, 729, 716, 715, 706, 709, 
  744, 751, 742, 737, 756, 755, 762, 765, 
  0, 8, 16, 24, 32, 40, 48, 56, 
  64, 72, 80, 88, 96, 104, 112, 120, 
  128, 136, 144, 152, 160, 168, 176, 184, 
  192, 200, 208, 216, 224, 232, 240, 248, 
  256, 264, 272, 280, 288, 296, 304, 312, 
  320, 328, 336, 344, 352, 360, 368, 376, 
  384, 392, 400, 408, 416, 424, 432, 440, 
  448, 456, 464, 472, 480, 488, 496, 504, 
  512, 520, 528, 536, 544, 552, 560, 568, 
  576, 584, 592, 600, 608, 616, 624, 632, 
  640, 648, 656, 664, 672, 680, 688, 696, 
  704, 712, 720, 728, 736, 744, 752, 760, 
  768, 776, 784, 792, 800, 808, 816, 824, 
  832, 840, 848, 856, 864, 872, 880, 888, 
  896, 904, 912, 920, 928, 936, 944, 952, 
  960, 968, 976, 984, 992, 1000, 1008, 1016, 
  1024, 1032, 1040, 1048, 1056, 1064, 1072, 1080, 
  1088, 1096, 1104, 1112, 1120, 1128, 1136, 1144, 
  1152, 1160, 1168, 1176, 1184, 1192, 1200, 1208, 
  1216, 1224, 1232, 1240, 1248, 1256, 1264, 1272, 
  1280, 1288, 1296, 1304, 1312, 1320, 1328, 1336, 
  1344, 1352, 1360, 1368, 1376, 1384, 1392, 1400, 
  1408, 1416, 1424, 1432, 1440, 1448, 1456, 1464, 
  1472, 1480, 1488, 1496, 1504, 1512, 1520, 1528, 
  1536, 1544, 1552, 1560, 1568, 1576, 1584, 1592, 
  1600, 1608, 1616, 1624, 1632, 1640, 1648, 1656, 
  1664, 1672, 1680, 1688, 1696, 1704, 1712, 1720, 
  1728, 1736, 1744, 1752, 1760, 1768, 1776, 1784, 
  1792, 1800, 1808, 1816, 1824, 1832, 1840, 1848, 
  1856, 1864, 1872, 1880, 1888, 1896, 1904, 1912, 
  1920, 1928, 1936, 1944, 1952, 1960, 1968, 1976, 
  1984, 1992, 2000, 2008, 2016, 2024, 2032, 2040, 
  0, 9, 18, 27, 36, 45, 54, 63, 
  72, 65, 90, 83, 108, 101, 126, 119, 
  144, 153, 130, 139, 180, 189, 166, 175, 
  216, 209, 202, 195, 252, 245, 238, 231, 
  288, 297, 306, 315, 260, 269, 278, 287, 
  360, 353, 378, 371, 332, 325, 350, 343, 
  432, 441, 418, 427, 404, 413, 390, 399, 
  504, 497, 490, 483, 476, 469, 462, 455, 
  576, 585, 594, 603, 612, 621, 630, 639, 
  520, 513, 538, 531, 556, 549, 574, 567, 
  720, 729, 706, 715, 756, 765, 742, 751, 
  664, 657, 650, 643, 700, 693, 686, 679, 
  864, 873, 882, 891, 836, 845, 854, 863, 
  808, 801, 826, 819, 780, 773, 798, 791, 
  1008, 1017, 994, 1003, 980, 989, 966, 975, 
  952, 945, 938, 931, 924, 917, 910, 903, 
  1152, 1161, 1170, 1179, 1188, 1197, 1206, 1215, 
  1224, 1217, 1242, 1235, 1260, 1253, 1278, 1271, 
  1040, 1049, 1026, 1035, 1076, 1085, 1062, 1071, 
  1112, 1105, 1098, 1091, 1148, 1141, 1134, 1127, 
  1440, 1449, 1458, 1467, 1412, 1421, 1430, 1439, 
  1512, 1505, 1530, 1523, 1484, 1477, 1502, 1495, 
  1328, 1337, 1314, 1323, 1300, 1309, 1286, 1295, 
  1400, 1393, 1386, 1379, 1372, 1365, 1358, 1351, 
  1728, 1737, 1746, 1755, 1764, 1773, 1782, 1791, 
  1672, 1665, 1690, 1683, 1708, 1701, 1726, 1719, 
  1616, 1625, 1602, 1611, 1652, 1661, 1638, 1647, 
  1560, 1553, 1546, 1539, 1596, 1589, 1582, 1575, 
  2016, 2025, 2034, 2043, 1988, 1997, 2006, 2015, 
  1960, 1953, 1978, 1971, 1932, 1925, 1950, 1943, 
  1904, 1913, 1890, 1899, 1876, 1885, 1862, 1871, 
  1848, 1841, 1834, 1827, 1820, 1813, 1806, 1799, 
  0, 10, 20, 30, 40, 34, 60, 54, 
  80, 90, 68, 78, 120, 114, 108, 102, 
  160, 170, 180, 190, 136, 130, 156, 150, 
  240, 250, 228, 238, 216, 210, 204, 198, 
  320, 330, 340, 350, 360, 354, 380, 374, 
  272, 282, 260, 270, 312, 306, 300, 294, 
  480, 490, 500, 510, 456, 450, 476, 470, 
  432, 442, 420, 430, 408, 402, 396, 390, 
  640, 650, 660, 670, 680, 674, 700, 694, 
  720, 730, 708, 718, 760, 754, 748, 742, 
  544, 554, 564, 574, 520, 514, 540, 534, 
  624, 634, 612, 622, 600, 594, 588, 582, 
  960, 970, 980, 990, 1000, 994, 1020, 1014, 
  912, 922, 900, 910, 952, 946, 940, 934, 
  864, 874, 884, 894, 840, 834, 860, 854, 
  816, 826, 804, 814, 792, 786, 780, 774, 
  1280, 1290, 1300, 1310, 1320, 1314, 1340, 1334, 
  1360, 1370, 1348, 1358, 1400, 1394, 1388, 1382, 
  1440, 1450, 1460, 1470, 1416, 1410, 1436, 1430, 
  1520, 1530, 1508, 1518, 1496, 1490, 1484, 1478, 
  1088, 1098, 1108, 1118, 1128, 1122, 1148, 1142, 
  1040, 1050, 1028, 1038, 1080, 1074, 1068, 1062, 
  1248, 1258, 1268, 1278, 1224, 1218, 1244, 1238, 
  1200, 1210, 1188, 1198, 1176, 1170, 1164, 1158, 
  1920, 1930, 1940, 1950, 1960, 1954, 1980, 1974, 
  2000, 2010, 1988, 1998, 2040, 2034, 2028, 2022, 
  1824, 1834, 1844, 1854, 1800, 1794, 1820, 1814, 
  1904, 1914, 1892, 1902, 1880, 1874, 1868, 1862, 
  1728, 1738, 1748, 1758, 1768, 1762, 1788, 1782, 
  1680, 1690, 1668, 1678, 1720, 1714, 1708, 1702, 
  1632, 1642, 1652, 1662, 1608, 1602, 1628, 1622, 
  1584, 1594, 1572, 1582, 1560, 1554, 1548, 1542, 
  0, 11, 22, 29, 44, 39, 58, 49, 
  88, 83, 78, 69, 116, 127, 98, 105, 
  176, 187, 166, 173, 156, 151, 138, 129, 
  232, 227, 254, 245, 196, 207, 210, 217, 
  352, 363, 374, 381, 332, 327, 346, 337, 
  312, 307, 302, 293, 276, 287, 258, 265, 
  464, 475, 454, 461, 508, 503, 490, 481, 
  392, 387, 414, 405, 420, 431, 434, 441, 
  704, 715, 726, 733, 748, 743, 762, 753, 
  664, 659, 654, 645, 692, 703, 674, 681, 
  624, 635, 614, 621, 604, 599, 586, 577, 
  552, 547, 574, 565, 516, 527, 530, 537, 
  928, 939, 950, 957, 908, 903, 922, 913, 
  1016, 1011, 1006, 997, 980, 991, 962, 969, 
  784, 795, 774, 781, 828, 823, 810, 801, 
  840, 835, 862, 853, 868, 879, 882, 889, 
  1408, 1419, 1430, 1437, 1452, 1447, 1466, 1457, 
  1496, 1491, 1486, 1477, 1524, 1535, 1506, 1513, 
  1328, 1339, 1318, 1325, 1308, 1303, 1290, 1281, 
  1384, 1379, 1406, 1397, 1348, 1359, 1362, 1369, 
  1248, 1259, 1270, 1277, 1228, 1223, 1242, 1233, 
  1208, 1203, 1198, 1189, 1172, 1183, 1154, 1161, 
  1104, 1115, 1094, 1101, 1148, 1143, 1130, 1121, 
  1032, 1027, 1054, 1045, 1060, 1071, 1074, 1081, 
  1856, 1867, 1878, 1885, 1900, 1895, 1914, 1905, 
  1816, 1811, 1806, 1797, 1844, 1855, 1826, 1833, 
  2032, 2043, 2022, 2029, 2012, 2007, 1994, 1985, 
  1960, 1955, 1982, 1973, 1924, 1935, 1938, 1945, 
  1568, 1579, 1590, 1597, 1548, 1543, 1562, 1553, 
  1656, 1651, 1646, 1637, 1620, 1631, 1602, 1609, 
  1680, 1691, 1670, 1677, 1724, 1719, 1706, 1697, 
  1736, 1731, 1758, 1749, 1764, 1775, 1778, 1785, 
  0, 12, 24, 20, 48, 60, 40, 36, 
  96, 108, 120, 116, 80, 92, 72, 68, 
  192, 204, 216, 212, 240, 252, 232, 228, 
  160, 172, 184, 180, 144, 156, 136, 132, 
  384, 396, 408, 404, 432, 444, 424, 420, 
  480, 492, 504, 500, 464, 476, 456, 452, 
  320, 332, 344, 340, 368, 380, 360, 356, 
  288, 300, 312, 308, 272, 284, 264, 260, 
  768, 780, 792, 788, 816, 828, 808, 804, 
  864, 876, 888, 884, 848, 860, 840, 836, 
  960, 972, 984, 980, 1008, 1020, 1000, 996, 
  928, 940, 952, 948, 912, 924, 904, 900, 
  640, 652, 664, 660, 688, 700, 680, 676, 
  736, 748, 760, 756, 720, 732, 712, 708, 
  576, 588, 600, 596, 624, 636, 616, 612, 
  544, 556, 568, 564, 528, 540, 520, 516, 
  1536, 1548, 1560, 1556, 1584, 1596, 1576, 1572, 
  1632, 1644, 1656, 1652, 1616, 1628, 1608, 1604, 
  1728, 1740, 1752, 1748, 1776, 1788, 1768, 1764, 
  1696, 1708, 1720, 1716, 1680, 1692, 1672, 1668, 
  1920, 1932, 1944, 1940, 1968, 1980, 1960, 1956, 
  2016, 2028, 2040, 2036, 2000, 2012, 1992, 1988, 
  1856, 1868, 1880, 1876, 1904, 1916, 1896, 1892, 
  1824, 1836, 1848, 1844, 1808, 1820, 1800, 1796, 
  1280, 1292, 1304, 1300, 1328, 1340, 1320, 1316, 
  1376, 1388, 1400, 1396, 1360, 1372, 1352, 1348, 
  1472, 1484, 1496, 1492, 1520, 1532, 1512, 1508, 
  1440, 1452, 1464, 1460, 1424, 1436, 1416, 1412, 
  1152, 1164, 1176, 1172, 1200, 1212, 1192, 1188, 
  1248, 1260, 1272, 1268, 1232, 1244, 1224, 1220, 
  1088, 1100, 1112, 1108, 1136, 1148, 1128, 1124, 
  1056, 1068, 1080, 1076, 1040, 1052, 1032, 1028, 
  0, 13, 26, 23, 52, 57, 46, 35, 
  104, 101, 114, 127, 92, 81, 70, 75, 
  208, 221, 202, 199, 228, 233, 254, 243, 
  184, 181, 162, 175, 140, 129, 150, 155, 
  416, 429, 442, 439, 404, 409, 398, 387, 
  456, 453, 466, 479, 508, 497, 486, 491, 
  368, 381, 362, 359, 324, 329, 350, 339, 
  280, 277, 258, 271, 300, 289, 310, 315, 
  832, 845, 858, 855, 884, 889, 878, 867, 
  808, 805, 818, 831, 796, 785, 774, 779, 
  912, 925, 906, 903, 932, 937, 958, 947, 
  1016, 1013, 994, 1007, 972, 961, 982, 987, 
  736, 749, 762, 759, 724, 729, 718, 707, 
  648, 645, 658, 671, 700, 689, 678, 683, 
  560, 573, 554, 551, 516, 521, 542, 531, 
  600, 597, 578, 591, 620, 609, 630, 635, 
  1664, 1677, 1690, 1687, 1716, 1721, 1710, 1699, 
  1768, 1765, 1778, 1791, 1756, 1745, 1734, 1739, 
  1616, 1629, 1610, 1607, 1636, 1641, 1662, 1651, 
  1592, 1589, 1570, 1583, 1548, 1537, 1558, 1563, 
  1824, 1837, 1850, 1847, 1812, 1817, 1806, 1795, 
  1864, 1861, 1874, 1887, 1916, 1905, 1894, 1899, 
  2032, 2045, 2026, 2023, 1988, 1993, 2014, 2003, 
  1944, 1941, 1922, 1935, 1964, 1953, 1974, 1979, 
  1472, 1485, 1498, 1495, 1524, 1529, 1518, 1507, 
  1448, 1445, 1458, 1471, 1436, 1425, 1414, 1419, 
  1296, 1309, 1290, 1287, 1316, 1321, 1342, 1331, 
  1400, 1397, 1378, 1391, 1356, 1345, 1366, 1371, 
  1120, 1133, 1146, 1143, 1108, 1113, 1102, 1091, 
  1032, 1029, 1042, 1055, 1084, 1073, 1062, 1067, 
  1200, 1213, 1194, 1191, 1156, 1161, 1182, 1171, 
  1240, 1237, 1218, 1231, 1260, 1249, 1270, 1275, 
  0, 14, 28, 18, 56, 54, 36, 42, 
  112, 126, 108, 98, 72, 70, 84, 90, 
  224, 238, 252, 242, 216, 214, 196, 202, 
  144, 158, 140, 130, 168, 166, 180, 186, 
  448, 462, 476, 466, 504, 502, 484, 490, 
  432, 446, 428, 418, 392, 390, 404, 410, 
  288, 302, 316, 306, 280, 278, 260, 266, 
  336, 350, 332, 322, 360, 358, 372, 378, 
  896, 910, 924, 914, 952, 950, 932, 938, 
  1008, 1022, 1004, 994, 968, 966, 980, 986, 
  864, 878, 892, 882, 856, 854, 836, 842, 
  784, 798, 780, 770, 808, 806, 820, 826, 
  576, 590, 604, 594, 632, 630, 612, 618, 
  560, 574, 556, 546, 520, 518, 532, 538, 
  672, 686, 700, 690, 664, 662, 644, 650, 
  720, 734, 716, 706, 744, 742, 756, 762, 
  1792, 1806, 1820, 1810, 1848, 1846, 1828, 1834, 
  1904, 1918, 1900, 1890, 1864, 1862, 1876, 1882, 
  2016, 2030, 2044, 2034, 2008, 2006, 1988, 1994, 
  1936, 1950, 1932, 1922, 1960, 1958, 1972, 1978, 
  1728, 1742, 1756, 1746, 1784, 1782, 1764, 1770, 
  1712, 1726, 1708, 1698, 1672, 1670, 1684, 1690, 
  1568, 1582, 1596, 1586, 1560, 1558, 1540, 1546, 
  1616, 1630, 1612, 1602, 1640, 1638, 1652, 1658, 
  1152, 1166, 1180, 1170, 1208, 1206, 1188, 1194, 
  1264, 1278, 1260, 1250, 1224, 1222, 1236, 1242, 
  1120, 1134, 1148, 1138, 1112, 1110, 1092, 1098, 
  1040, 1054, 1036, 1026, 1064, 1062, 1076, 1082, 
  1344, 1358, 1372, 1362, 1400, 1398, 1380, 1386, 
  1328, 1342, 1324, 1314, 1288, 1286, 1300, 1306, 
  1440, 1454, 1468, 1458, 1432, 1430, 1412, 1418, 
  1488, 1502, 1484, 1474, 1512, 1510, 1524, 1530, 
  0, 15, 30, 17, 60, 51, 34, 45, 
  120, 119, 102, 105, 68, 75, 90, 85, 
  240, 255, 238, 225, 204, 195, 210, 221, 
  136, 135, 150, 153, 180, 187, 170, 165, 
  480, 495, 510, 497, 476, 467, 450, 461, 
  408, 407, 390, 393, 420, 427, 442, 437, 
  272, 287, 270, 257, 300, 291, 306, 317, 
  360, 359, 374, 377, 340, 347, 330, 325, 
  960, 975, 990, 977, 1020, 1011, 994, 1005, 
  952, 951, 934, 937, 900, 907, 922, 917, 
  816, 831, 814, 801, 780, 771, 786, 797, 
  840, 839, 854, 857, 884, 891, 874, 869, 
  544, 559, 574, 561, 540, 531, 514, 525, 
  600, 599, 582, 585, 612, 619, 634, 629, 
  720, 735, 718, 705, 748, 739, 754, 765, 
  680, 679, 694, 697, 660, 667, 650, 645, 
  1920, 1935, 1950, 1937, 1980, 1971, 1954, 1965, 
  2040, 2039, 2022, 2025, 1988, 1995, 2010, 2005, 
  1904, 1919, 1902, 1889, 1868, 1859, 1874, 1885, 
  1800, 1799, 1814, 1817, 1844, 1851, 1834, 1829, 
  1632, 1647, 1662, 1649, 1628, 1619, 1602, 1613, 
  1560, 1559, 1542, 1545, 1572, 1579, 1594, 1589, 
  1680, 1695, 1678, 1665, 1708, 1699, 1714, 1725, 
  1768, 1767, 1782, 1785, 1748, 1755, 1738, 1733, 
  1088, 1103, 1118, 1105, 1148, 1139, 1122, 1133, 
  1080, 1079, 1062, 1065, 1028, 1035, 1050, 1045, 
  1200, 1215, 1198, 1185, 1164, 1155, 1170, 1181, 
  1224, 1223, 1238, 1241, 1268, 1275, 1258, 1253, 
  1440, 1455, 1470, 1457, 1436, 1427, 1410, 1421, 
  1496, 1495, 1478, 1481, 1508, 1515, 1530, 1525, 
  1360, 1375, 1358, 1345, 1388, 1379, 1394, 1405, 
  1320, 1319, 1334, 1337, 1300, 1307, 1290, 1285, 
];

// 

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn mull_as_shift() {
	// multiplying by a power of x is equivalent to left-shift
	for byte in 0u8..=255 {
	    for bit in 0..4 {
		assert_eq!(rmull(byte, 1 << bit),
			   (byte as u16) << bit);
	    }
	}
    }

    #[test]
    fn lmull_vs_rmull() {
	// lmull values should simply be rmull values << 4
	for byte in 0u8..=255 {
	    for nibble in 0u8..16 {
		assert_eq!(lmull(byte, nibble),
			   rmull(byte, nibble) << 4);
	    }
	}
    }

    #[test]
    fn lrmull_vs_lmull_rmull() {
	for byte in 0u8..=255 {
	    for l in 0u8..16 {
		for r in 0u8..16 {
		    let combined : u8 = (l << 4) | r;
		    assert_eq!(lrmull(byte, combined),
			       lmull(byte, l) ^ rmull(byte, r));
		}
	    }
	}	
    }
}
