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
	union {
		real s[3];
		struct { real s0, s1, s2; };
		struct { real x, y, z; };
	};
};

struct real4 {
	union {
		real s[4];
		struct { real s0, s1, s2, s3; };
	//	struct { real t, x, y, z; };
	};
};

struct real3s3 {
	union {
		real s[6];
		struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
	//	struct { real xx, xy, xz, yy, yz, zz; };
	};
};

struct real4s4 {
	union {
		real s[10];
		struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
	//	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
	};
};

struct real4x4 {
	union {
		real4 s[4];
		struct { real4 s0, s1, s2, s3; };
	//	struct { real4 t, x, y, z; };
	};
};

struct real4x4s4 {
	union {
		real4s4 s[4];
		struct { real4s4 s0, s1, s2, s3; };
	//	struct { real4s4 t, x, y, z; };
	};
};

struct real4x4x4 {
	union {
		real4x4 s[4];
		struct { real4x4 s0, s1, s2, s3; };
	//	struct { real4x4 t, x, y, z; };
	};
};

struct real4x4x4s4 {
	union {
		real4x4s4 s[4];
		struct { real4x4s4 s0, s1, s2, s3; };
	//	struct { real4x4s4 t, x, y, z; };
	};
};

struct real4s4x4s4 {
	union {
		real4s4 s[10];
		struct { real4s4 s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };
	//	struct { real4s4 tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
	};
};

struct real4x4x4x4 {
	union {
		real4x4x4 s[4];
		struct { real4x4x4 s0, s1, s2, s3; };
	//	struct { real4x4x4 t, x, y, z; };
	};
};

// one option: I could define structs using luajit ffi struct lib
// then I can auto validate their alignment with opencl
// but I'd have to hack in ctors etc in efe.lua
// another option: put them here, and put equivalent code in luajit ffi only
// I think I'll do that, so I can define these as proper C++ classes...

#endif


