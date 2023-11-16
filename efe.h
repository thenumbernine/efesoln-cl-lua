typedef union {
	real s[3];
	struct { real s0, s1, s2; };
	struct { real x, y, z; };
} real3;

typedef union {
	real s[6];
	struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
	struct { real xx, xy, xz, yy, yz, zz; };
} real3s3;

typedef union {
	real s[10];
	struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
} real4s4;

typedef union {
	real s[16];
	struct { real4 s0, s1, s2, s3; };
	struct { real4 t, x, y, z; };
} real4x4;

typedef union {
	real4s4 s[4];
	struct { real4s4 s0, s1, s2, s3; };
	struct { real4s4 t, x, y, z; };
} real4x4s4;

typedef union {
	real4x4 s[4];
	struct { real4x4 s0, s1, s2, s3; };
	struct { real4x4 t, x, y, z; };
} real4x4x4;

typedef union {
	real4x4s4 s[4];
	struct { real4x4s4 s0, s1, s2, s3; };
	struct { real4x4s4 t, x, y, z; };
} real4x4x4s4;

typedef union {
	real4s4 s[10];
	struct { real4s4 s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };
	struct { real4s4 tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
} real4s4x4s4;
