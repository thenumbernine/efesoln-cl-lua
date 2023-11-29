#pragma once

// this is the header that clcpp files will use
// it has to be included after autogen.h which defines real
// but the #include can't go in here cuz luajit ffi doesn't like them

// stuck on:
// ../../cpp/Common/include/Common/Variadic.h:6:10: fatal error: 'cstddef' file not found
// how to get cstddef into clcpp ...
#if 1

#include "Tensor/Vector.h"

using real3 = Tensor::vec<real, 3>;
using real4 = Tensor::vec<real, 4>;
using real3x3 = Tensor::mat<real, 3, 3>;
using real3s3 = Tensor::sym<real, 3>;
using real4s4 = Tensor::sym<real, 4>;
using real4x4 = Tensor::mat<real, 4, 4>;
using real4x4s4 = Tensor::tensorx<real, 4, -'s', 4>;
using real4x4x4 = Tensor::tensorr<real, 4, 3>;
using real4x4x4s4 = Tensor::tensorx<real, 4, 4, -'s', 4>;
using real4s4x4s4 = Tensor::tensorx<real, -'s', 4, -'s', 4>;
using real4x4x4x4 = Tensor::tensorr<real, 4, 4>;

#endif

#if 0
struct real3 {
	union {
		real s[3] = { real(0) };
		struct { real s0, s1, s2; };
		struct { real x, y, z; };
	};

	real3() {}

	real3(
		real const s0_,
		real const s1_,
		real const s2_
	) : s0(s0_),
		s1(s1_),
		s2(s2_)
	{}
};

struct real4 {
	union {
		real s[4] = { real(0) };
		struct { real s0, s1, s2, s3; };
	//	struct { real t, x, y, z; };
	};

	real4() {}

	real4(
		real const s0_,
		real const s1_,
		real const s2_,
		real const s3_
	) : s0(s0_),
		s1(s1_),
		s2(s2_),
		s3(s3_)
	{}
};

struct real3s3 {
	union {
		real s[6] = { real(0) };
		struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
	//	struct { real xx, xy, xz, yy, yz, zz; };
	};

	//notice this is a different order than Tensor, which appends new dimensions to the end of the dense triangular matrix, not to the beginning like this
	real3s3() {}

	real3s3(
		real const s00_,
		real const s01_,
		real const s02_,
		real const s11_,
		real const s12_,
		real const s22_
	) : s00(s00_),
		s01(s01_),
		s02(s02_),
		s11(s11_),
		s12(s12_),
		s22(s22_)
	{}
};

struct real4s4 {
	union {
		real s[10] = { real(0) };
		struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
	//	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
	};
	
	real4s4() {}

	real4s4(
		real const s00_,
		real const s01_,
		real const s02_,
		real const s03_,
		real const s11_,
		real const s12_,
		real const s13_,
		real const s22_,
		real const s23_,
		real const s33_
	) : s00(s00_),
		s01(s01_),
		s02(s02_),
		s03(s03_),
		s11(s11_),
		s12(s12_),
		s13(s13_),
		s22(s22_),
		s23(s23_),
		s33(s33_)
	{}
};

struct real4x4 {
	union {
		real4 s[4] = { real4() };
		struct { real4 s0, s1, s2, s3; };
	//	struct { real4 t, x, y, z; };
	};

	real4x4() {}

	real4x4(
		real4 const s0_,
		real4 const s1_,
		real4 const s2_,
		real4 const s3_
	) : s0(s0_),
		s1(s1_),
		s2(s2_),
		s3(s3_)
	{}
};

struct real4x4s4 {
	union {
		real4s4 s[4] = { real4s4() };
		struct { real4s4 s0, s1, s2, s3; };
	//	struct { real4s4 t, x, y, z; };
	};

	real4x4s4() {}

	real4x4s4(
		real4s4 const s0_,
		real4s4 const s1_,
		real4s4 const s2_,
		real4s4 const s3_
	) : s0(s0_),
		s1(s1_),
		s2(s2_),
		s3(s3_)
	{}
};

struct real4x4x4 {
	union {
		real4x4 s[4] = { real4x4() };
		struct { real4x4 s0, s1, s2, s3; };
	//	struct { real4x4 t, x, y, z; };
	};

	real4x4x4() {}

	real4x4x4(
		real4x4 const s0_,
		real4x4 const s1_,
		real4x4 const s2_,
		real4x4 const s3_
	) : s0(s0_),
		s1(s1_),
		s2(s2_),
		s3(s3_)
	{}
};

struct real4x4x4s4 {
	union {
		real4x4s4 s[4] = { real4x4s4() };
		struct { real4x4s4 s0, s1, s2, s3; };
	//	struct { real4x4s4 t, x, y, z; };
	};

	real4x4x4s4() {}

	real4x4x4s4(
		real4x4s4 const s0_,
		real4x4s4 const s1_,
		real4x4s4 const s2_,
		real4x4s4 const s3_
	) : s0(s0_),
		s1(s1_),
		s2(s2_),
		s3(s3_)
	{}
};

struct real4s4x4s4 {
	union {
		real4s4 s[10] = { real4s4() };
		struct { real4s4 s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };
	//	struct { real4s4 tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
	};

	real4s4x4s4() {}

	real4s4x4s4(
		real4s4 const s00_,
		real4s4 const s01_,
		real4s4 const s02_,
		real4s4 const s03_,
		real4s4 const s11_,
		real4s4 const s12_,
		real4s4 const s13_,
		real4s4 const s22_,
		real4s4 const s23_,
		real4s4 const s33_
	) : s00(s00_),
		s01(s01_),
		s02(s02_),
		s03(s03_),
		s11(s11_),
		s12(s12_),
		s13(s13_),
		s22(s22_),
		s23(s23_),
		s33(s33_)
	{}
};

struct real4x4x4x4 {
	union {
		real4x4x4 s[4] = { real4x4x4() };
		struct { real4x4x4 s0, s1, s2, s3; };
	//	struct { real4x4x4 t, x, y, z; };
	};
	
	real4x4x4x4() {}

	real4x4x4x4(
		real4x4x4 const s0_,
		real4x4x4 const s1_,
		real4x4x4 const s2_,
		real4x4x4 const s3_
	) : s0(s0_),
		s1(s1_),
		s2(s2_),
		s3(s3_)
	{}
};

// one option: I could define structs using luajit ffi struct lib
// then I can auto validate their alignment with opencl
// but I'd have to hack in ctors etc in efe.lua
// another option: put them here, and put equivalent code in luajit ffi only
// I think I'll do that, so I can define these as proper C++ classes...

#endif
