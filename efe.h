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

typedef struct {
	real alpha;
	real3 betaU;
	sym3 gammaLL;
} gPrim_t;

typedef struct {
	//source terms:
<? if solver.body.useMatter then ?>
	real rho;	//matter density
	real P;		//pressure ... due to matter.  TODO what about magnetic pressure?
	real eInt;	//specific internal energy
<?	if solver.body.useVel then ?>
	real3 v;	//3-vel (upper, spatial)
<? 	end
end ?>
<? if solver.body.useEM then ?>
	<? if solver.useFourPotential then ?>
	real4 JL;	//4-current: rho, j  use to solve ...
	real4 AL;	//4-potential: phi, A
	<? else ?>
//this needs to be lienar solved for ... but it's an easy problem (at least when geometry is flat)
//	real chargeDensity;
//	TensorUsub currentDensity;	//TODO how does this relate to matter density?

//in the mean time ...
	real3 E, B;	//upper, spatial
	<? end ?>
<? end ?>
} <?=TPrim_t?>;
