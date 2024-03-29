// this is the luajit version of the structures
// they have to have the same fields
// then i autogen code to make sure the field offsets and sizes match

typedef union {
	real s[3];
	struct { real s0, s1, s2; };
	struct { real x, y, z; };
} real3;

typedef union {
	real s[4];
	struct { real s0, s1, s2, s3; };
//	struct { real t, x, y, z; };
} real4;

typedef union {
	real s[6];
	struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
//	struct { real xx, xy, xz, yy, yz, zz; };
} real3s3;

typedef union {
	real s[10];
	struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
//	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
} real4s4;

typedef union {
	real4 s[4];
	struct { real4 s0, s1, s2, s3; };
//	struct { real4 t, x, y, z; };
} real4x4;

typedef union {
	real4s4 s[4];
	struct { real4s4 s0, s1, s2, s3; };
//	struct { real4s4 t, x, y, z; };
} real4x4s4;

typedef union {
	real4x4 s[4];
	struct { real4x4 s0, s1, s2, s3; };
//	struct { real4x4 t, x, y, z; };
} real4x4x4;

typedef union {
	real4x4s4 s[4];
	struct { real4x4s4 s0, s1, s2, s3; };
//	struct { real4x4s4 t, x, y, z; };
} real4x4x4s4;

typedef union {
	real4s4 s[10];
	struct { real4s4 s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };
//	struct { real4s4 tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
} real4s4x4s4;

typedef union {
	real4x4x4 s[4];
	struct { real4x4x4 s0, s1, s2, s3; };
//	struct { real4x4x4 t, x, y, z; };
} real4x4x4x4;


// one option: I could define structs using luajit ffi struct lib
// then I can auto validate their alignment with opencl
// but I'd have to hack in ctors etc in efe.lua
// another option: put them here, and put equivalent code in luajit ffi only
// I think I'll do that, so I can define these as proper C++ classes...
