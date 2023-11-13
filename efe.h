typedef union {
	real s[3];
	struct { real s0, s1, s2; };
	struct { real x, y, z; };
} real3;

typedef union {
	real s[6];
	struct { real s00, s01, s02, s11, s12, s22; };	//useful for templated code
	struct { real xx, xy, xz, yy, yz, zz; };
} sym3;

typedef union {
	real s[10];
	struct { real s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };	//useful for templated code
	struct { real tt, tx, ty, tz, xx, xy, xz, yy, yz, zz; };
} sym4;

typedef union {
	sym4 s[4];
	struct { sym4 s0, s1, s2, s3; };
} tensor_4sym4;

typedef union {
	tensor_4sym4 s[4];
	struct { tensor_4sym4 s0, s1, s2, s3; };
} tensor_44sym4;

typedef union {
	sym4 s[10];
	struct { sym4 s00, s01, s02, s03, s11, s12, s13, s22, s23, s33; };
} tensor_sym4sym4;
