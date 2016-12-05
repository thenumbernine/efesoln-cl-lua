<? if real == 'double' then ?>
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
<? end ?>

typedef <?=real?> real;
typedef <?=real?>2 real2;
typedef <?=real?>4 real4;

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
	real s[4*4];
	struct { real4 s0, s1, s2, s3; };
	struct {
<? for a=0,3 do
	for b=0,3 do ?>
		real s<?=a?><?=b?>;
<?	end
end ?>
	};
} mat4;

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

typedef struct {
	real alpha;
	real3 betaU;
	sym3 gammaLL;
} gPrim_t;

typedef struct {
	//source terms:
	real rho;	//matter density
	real P;		//pressure ... due to matter.  TODO what about magnetic pressure?
	real eInt;	//specific internal energy

	real3 v;	//3-vel (upper, spatial)

//this needs to be lienar solved for ... but it's an easy problem (at least when geometry is flat)
//	real chargeDensity;
//	TensorUsub currentDensity;	//TODO how does this relate to matter density?

//in the mean time ...
	real3 E, B;	//upper, spatial

} TPrim_t;
