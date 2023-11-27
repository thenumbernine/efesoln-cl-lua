#pragma once

// this is the header that clcpp files will use
// it has to be included after autogen.h which defines real
// but the #include can't go in here cuz luajit ffi doesn't like them

// stuck on:
// ../../cpp/Common/include/Common/Variadic.h:6:10: fatal error: 'cstddef' file not found
// how to get cstddef into clcpp ...
#if 0
#include "Tensor/Vector.h"
using real3 = Tensor::vec<real, 3>;
using real4 = Tensor::vec<real, 4>;
#endif

#if 1

struct real3 {
	//union {
		real s[3];
	//	struct { real s0, s1, s2; };
	//	struct { real x, y, z; };
	//};

	real3() {
		for (int i = 0; i < 3; ++i) {
			s[i] = real(0);
		}
	}

	real3(
		real const s0_,
		real const s1_,
		real const s2_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
	}
};

struct real4 {
	//union {
	real s[4];
	//	struct { real s0, s1, s2, s3; };
	//	struct { real t, x, y, z; };
	//};

	real4() {
		for (int i = 0; i < 4; ++i) {
			s[i] = real(0);
		}
	}

	real4(
		real const s0_,
		real const s1_,
		real const s2_,
		real const s3_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
	}
};

struct real3s3 {
	//union {
	real s[6];
	//	struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
	//	struct { real xx, xy, xz, yy, yz, zz; };
	//};

	real3s3() {
		for (int i = 0; i < 6; ++i) {
			s[i] = real(0);
		}
	}

	real3s3(
		real const s0_,
		real const s1_,
		real const s2_,
		real const s3_,
		real const s4_,
		real const s5_
	) { 
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
		s[4] = s4_;
		s[5] = s5_;
	}
};

struct real4s4 {
	//union {
	real s[10];
	//	struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
	//	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
	//};
	
	real4s4() {
		for (int i = 0; i < 10; ++i) {
			s[i] = real(0);
		}
	}

	real4s4(
		real const s0_,
		real const s1_,
		real const s2_,
		real const s3_,
		real const s4_,
		real const s5_,
		real const s6_,
		real const s7_,
		real const s8_,
		real const s9_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
		s[4] = s4_;
		s[5] = s5_;
		s[6] = s6_;
		s[7] = s7_;
		s[8] = s8_;
		s[9] = s9_;
	}
};

struct real4x4 {
	//union {
	real4 s[4];
	//	struct { real4 s0, s1, s2, s3; };
	//	struct { real4 t, x, y, z; };
	//};

	real4x4() {
		for (int i = 0; i < 4; ++i) {
			s[i] = real4();
		}
	}

	real4x4(
		real4 const & s0_,
		real4 const & s1_,
		real4 const & s2_,
		real4 const & s3_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
	}
};

struct real4x4s4 {
	//union {
	real4s4 s[4];
	//	struct { real4s4 s0, s1, s2, s3; };
	//	struct { real4s4 t, x, y, z; };
	//};

	real4x4s4() {
		for (int i = 0; i < 4; ++i) {
			s[i] = real4s4();
		}
	}

	real4x4s4(
		real4s4 const & s0_,
		real4s4 const & s1_,
		real4s4 const & s2_,
		real4s4 const & s3_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
	}
};

struct real4x4x4 {
	//union {
	real4x4 s[4];
	//	struct { real4x4 s0, s1, s2, s3; };
	//	struct { real4x4 t, x, y, z; };
	//};

	real4x4x4() {
		for (int i = 0; i < 4; ++i) {
			s[i] = real4x4();
		}
	}

	real4x4x4(
		real4x4 const & s0_,
		real4x4 const & s1_,
		real4x4 const & s2_,
		real4x4 const & s3_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
	}
};

struct real4x4x4s4 {
	//union {
	real4x4s4 s[4];
	//	struct { real4x4s4 s0, s1, s2, s3; };
	//	struct { real4x4s4 t, x, y, z; };
	//};

	real4x4x4s4() {
		for (int i = 0; i < 4; ++i) {
			s[i] = real4x4s4();
		}
	}

	real4x4x4s4(
		real4x4s4 const & s0_,
		real4x4s4 const & s1_,
		real4x4s4 const & s2_,
		real4x4s4 const & s3_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
	}
};

struct real4s4x4s4 {
	//union {
	real4s4 s[10];
	//	struct { real4s4 s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };
	//	struct { real4s4 tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
	//};

	real4s4x4s4() {
		for (int i = 0; i < 10; ++i) {
			s[i] = real4s4();
		}
	}

	real4s4x4s4(
		real4s4 const & s0_,
		real4s4 const & s1_,
		real4s4 const & s2_,
		real4s4 const & s3_,
		real4s4 const & s4_,
		real4s4 const & s5_,
		real4s4 const & s6_,
		real4s4 const & s7_,
		real4s4 const & s8_,
		real4s4 const & s9_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
		s[4] = s4_;
		s[5] = s5_;
		s[6] = s6_;
		s[7] = s7_;
		s[8] = s8_;
		s[9] = s9_;
	}
};

struct real4x4x4x4 {
	//union {
	real4x4x4 s[4];
	//	struct { real4x4x4 s0, s1, s2, s3; };
	//	struct { real4x4x4 t, x, y, z; };
	//};
	
	real4x4x4x4() {
		for (int i = 0; i < 4; ++i) {
			s[i] = real4x4x4();
		}
	}

	real4x4x4x4(
		real4x4x4 const & s0_,
		real4x4x4 const & s1_,
		real4x4x4 const & s2_,
		real4x4x4 const & s3_
	) {
		s[0] = s0_;
		s[1] = s1_;
		s[2] = s2_;
		s[3] = s3_;
	}
};

// one option: I could define structs using luajit ffi struct lib
// then I can auto validate their alignment with opencl
// but I'd have to hack in ctors etc in efe.lua
// another option: put them here, and put equivalent code in luajit ffi only
// I think I'll do that, so I can define these as proper C++ classes...

#endif
